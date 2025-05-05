%%INIT
clear all; 
close all; 
clc; 
 
[x fs]= audioread('Drum+Bass.wav'); 
% 
% x=s; 
%% STFT 
FFTSIZE = 1024; 
HOPSIZE = 256; 
WINDOWSIZE = 512; 
 
X = myspectrogram(x,FFTSIZE,fs,hann(WINDOWSIZE),-HOPSIZE); 
V = abs(X(1:(FFTSIZE/2+1),:)); %Non-negative 
F = size(V,1); 
T = size(V,2); 
 
%sound(x,fs) 
%% NMF 
 
K = 3; %number of basis vectors 
MAXITER = 100; %total number of iterations to run 
 
%% TODO:NMF 
 
[W,H] = nmf(V,K,MAXITER); 
 
%% ISTFT / RECONSTRUCTION METHOD 1(SYNTHESIS) 
 
% get the mixture phase 
 
phi = angle(X); 
 
for i=1:K 
    XmagHat = W(:,i)*H(i,:); 
    % create upper half of frequency before istft 
    XmagHat = [XmagHat; conj(XmagHat(end-1:-1:2,:))]; 
    % Multiply with phase 
    XHat = XmagHat.*exp(1i*phi); 
    xhat1(:,i) = real(invmyspectrogram(XHat,HOPSIZE))'; 
end 
% sound(xhat1(:,1),fs) 
% sound(xhat1(:,2),fs) 
% sound(xhat1(:,3),fs) 
%% ISTFT / RECONSTRUCTION METHOD2 (FILTERING) 
% get the mixture phase 
phi = angle(X); 
for i=1:K 
    % create masking filture 
    Mask = (W(:,i)*H(i,:)./(W*H)); 
    % filter 
    XmagHat = V.*Mask; 
    %create upper half of frequency before istft 
    XmagHat = [XmagHat;conj(XmagHat(end-1:-1:2,:))]; 
    % Multiply with phase 
    XHat = XmagHat.*exp(1i*phi); 
    % create upper half of frequency before istft 
    xhat2(:,i) = real(invmyspectrogram(XmagHat.*exp(1i*phi),HOPSIZE))'; 
end 
% sound(xhat1(:,1),fs) 
% sound(xhat1(:,2),fs) 
% sound(xhat1(:,3),fs)

function [W,H] = nmf(V,K,MAXITER)
F = size(V,1);
T = size(V,2);
rand('seed',0)
W = 1+rand(F,K);
% W = W./repmat(sum(W),F,1);
H = 1+rand(K,T);
ONES = ones(F,T);
for i = 1:MAXITER    
    % update activations
    H = H.*(W'*(V./(W*H+eps)))./(W'*ONES);    
    % update dictionaries
    W = W.*((V./(W*H+eps))*H')./(ONES*H');
end
% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;


function X = myspectrogram(x,nfft,fs,window,noverlap,doplot,dbdown)
if nargin<7,dbdown=100;end
if nargin<6 doplot=0;end
if nargin<5 noverlap=256;end
if nargin<3,fs=1;end
if nargin<2,nfft=2048;end
x=x(:); %make sure it's a column
M=length(window);
if(M<2)error(...
        'myspectogram: Expect complete window, now just its length');
end;
if(M<2)error(...
        'myspectogram: Expect complete window, now just its length');
end;
if length(x)<M%zero-pad to fill a window:
    x=[x;zeros(M-length(x),1)];
end
Modd = mod(M,2);%0 if M even ,1 if odd
Mo2=(M-Modd)/2;
w=window(:);%Make sure it's a column
if noverlap<0
    nhop=-noverlap;
    noverlap = M-nhop;
else
    nhop = M-noverlap;
end
nx = length(x);
%nframes = 1+floor((nx-noverlap)/nhop);
nframes = 1+ceil(nx/nhop);
X = zeros(nfft,nframes);%allocate output spectogram
zp = zeros(nfft-M,1);%zero-paddingfor each FFT
xframe = zeros(M,1);
xoff=0-Mo2;%input time offset = half a frame
for m=1:nframes
    %M,Mo2,xoff,nhop
    if xoff<0
        xframe(1:xoff+M)=x(1:xoff+M);%partial input data frame
    else
        if xoff+M>nx
            xframe=[x(xoff+1:nx);zeros(xoff+M-nx,1)];
        else
            xframe=x(xoff+1:xoff+M);%input data frame
        end
    end
    xw = w.*xframe,%Apply window
    xwzp = [xw(Mo2+1:M);zp;xw(1:Mo2)];
    X(:,m)=fft(xwzp);
    xoff=xoff + nhop; % advancein[ut offset by hopsize
end
if(nargout==0)|doplot
    t=(0:nframes-1)*nhop/fs;
    f=0.001*(0:nfft-1)*fs/nfft;
    xdb = 20*log10(abs(X));
    Xmax = max(max(Xdb));
    %clip lower limit to -dbdown dB so nulls don't dominate:
    clipvals=[Xmax-dbdown,Xmax];
    imagesct(t,fXdb,clipvals);
    %grid;
    axis('xy');
    colormap(jet);
    xlabel('Time(sec)');
    ylabel('Freq(kHz)');
end


function a = invmyspectrogram(b,hop)
[nfft,nframes] = size(b);
No2 = nfft/2; % nfft assumed even
a = zeros(1, nfft+(nframes-1)*hop);
xoff = 0-No2; % output time offset = half of FFT size
for col = 1:nframes
    fftframe = b(:,col);    
    xzp = ifft(fftframe);
    % xzp = real(xzp); % if signal known to be real    
    x = [xzp(nfft-No2+1:nfft);xzp(1:No2)];
    if xoff<0 %FFT's "negative-time indices" are out of range        
        ix = 1:xoff+nfft;
        a(ix) = a(ix) + x(1-xoff:nfft)'; % partial frames out    
    else
        ix = xoff + 1:xoff+nfft;        
        a(ix) = a(ix) + x'; % overlap-add reconstruction
    end
    xoff = xoff + hop;
end

end
end
end

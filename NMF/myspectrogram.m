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

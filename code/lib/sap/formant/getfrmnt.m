function spectData = getfrmnt(frame, fs, lpcOrder, frmntmethod, plotopt)

% GETFRMNT Get Formant information : frequency and bandwidth
%    Usage : spectData = getfrmnt(frame, fs, lpcOrder, frmntmethod, plotopt)
%    The spectData is a structural variable inclusive several fields.
%    These fields are frmntfrq, frmntbw, pitch, fftspect, lpcspect, A, presig, poles, etc.
%    The frmntmethod presents several well known approaches to implement the procedure.
%    The lpcOrder usaully is set 12.
%
%    Cheng-Yuan Lin, 2003, January, 20.

if nargin < 5, plotopt = 0; end;
if nargin < 4, frmntmethod = 2; end;
if nargin < 3, lpcOrder = 12; end;
if nargin==0; selfdemo; return; end;

% To get spectrum information and related LPC poles solutions.
[fftspect, lpcspect, freqscale, A, presig, poles] = getspect(frame, fs, lpcOrder);
spectfactor = (fs/2)/length(fftspect);
Hzfactor    = fs/(2*pi);

switch frmntmethod,
case 1, %方法1 : W = arg(Z), B = -2*log(|Z|),  
   W = Hzfactor*angle(poles);
   B = Hzfactor*-2*log(abs(poles)); 
case 2, %方法2 : W 由 LPC spectrum 的 local maximum求取, B = 50 + 0.005*W
   [maxIndex, value] = localmax(lpcspect);
   W = maxIndex*spectfactor;
   B = 50 + 0.005*W;
end;

%由上面所求得的W, B再經一些formant基本定義的rules判斷是否合法
%index = find(W>90 & W<fs/2 & B<600);
index = find(W>90 & W<fs/2);
W = W(index);
B = B(index);
[W, ind] = sort(W);
B = B(ind);

%只保留前三個formant: F1, F2, F3
frmntfrq = W(1: min(3, length(W)));
frmntbw = B(1: min(3, length(B)));

%Computing pitch
pp = setPP;
pp.frameSize = length(frame);
pitch =  ptrackfcn(frame, fs, 0, pp);

%Assignment of spectData.
spectData.frmntfrq = frmntfrq;
spectData.frmntbw = frmntbw;
spectData.freqscale = freqscale;
spectData.pitch = pitch;
spectData.fftspect = fftspect;
spectData.lpcspect = lpcspect;
spectData.A = A;
spectData.presig = presig;
spectData.poles = poles;


%畫圖
if plotopt,
   figure('Position', [50 50 760 620]);
   subplot(2,1,1);
   plot(frame);
   title('Waveform');
   xlabel('Sample Points');
   ylabel('Magnitude');
   subplot(2,1,2);
   plot(freqscale, fftspect, 'b.-');
   minfft = min([fftspect lpcspect]) - 5;
   maxfft = max([fftspect lpcspect]) + 5;
   hold on;
   plot(freqscale, lpcspect, 'g+-');
   A = line([pitch pitch], [minfft maxfft]);
   set(A, 'linewidth', 2, 'color', 'black', 'linestyle','--');
   for k = 1 : length(frmntfrq),
      A = line([frmntfrq(k) frmntfrq(k)], [minfft maxfft]);
      set(A, 'linewidth', 2, 'color', 'red');
   end;
   legend('Fourier Spectrum', 'LPC Spectrum', 'F0 (pitch)', 'Formant Position');
   title('Spectrum');
   xlabel('Frequency (Hz)');
   ylabel('Log Magnitude (dB)');
   axis tight;
end;


function selfdemo
[y, fs] = wavread('test.wav');
y = y - mean(y);
%sound(y, fs);
lpcOrder = 12;
frmntmethod = 1; %best
plotopt = 1;
if fs>800,
   frameSize = 512;
else
   frameSize = 256;
end;
overlap   = (frameSize/2);
frmdY = buffer2(y, frameSize, overlap);
for k = 1 : size(frmdY, 2),
   frame = frmdY(:, k).*hamming(frameSize);
   spectData = getfrmnt(frame, fs, lpcOrder, frmntmethod, plotopt);
   pause;
   close all;
end;

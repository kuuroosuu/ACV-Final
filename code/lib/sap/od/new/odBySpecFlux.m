function [onsetTime, onsetFrameIdx, onsetSampleIdx] = odBySpecFlux(y,fs,nbits,odParam,plotOpt,groundTruth)
% Perform onset detection (OD) for a given music signal
% odBySpecFlux: OD based on Spectral Flux only 
%   Usage onset = odBySpecFlux(y,fs,odParam,plotOpt)
%       onset: detected onset time information (unit: sec)
%       y: input music signals
%       fs: sampling rate
%       odParam: parameters for onset detection
%       plotOpt: 0 for silence operation, 1 for plotting
%   
%   Example:
%       waveFile = 'little_star_piano2.wav';
%       [y,fs,nbits] = wavread(waveFile);
%		odParam=odParamSet2(fs);
%		plotOpt=1;
%		onsetTime=odBySpecFlux(y, fs, nbits, odParam, plotOpt);
% 
%   Jia-Min Ren, 20091104, 20091124

if nargin < 1, selfdemo, return; end
if nargin < 2, fs = 16000; end,
if nargin < 3, nbits=16; end
if nargin < 4, odParam = odParamSet2(fs); end,
if nargin < 5, plotOpt=0; end,

frameSize = odParam.frameSize;   % frame size
overlap = odParam.overlap;       % overlap
nfft = odParam.nfft;             % window size (resolution)

% y = y-mean(y);  % DC removal

frameMat=buffer2(y, frameSize, overlap);
frameNum=size(frameMat, 2);
SF=zeros(1, frameNum);       % detection function (spectral flux)

% DC removal
for i=1:frameNum,
    frameMat(:,i) = frameMat(:,i)-mean(frameMat(:,1));
end

% step 1. pre-processing 
% (multi-band processing, transient/steady-state separation ect...)
% (optional, not finished yet)

% step 2. reduction
% (using a single parameter to represent a frame)
% perform hamming windowing
hammWin = repmat(hamming(frameSize), 1, frameNum);
frameMat = frameMat .* hammWin;
% then calculate spectral flux from each frequency bin
magSpecMat = abs(fft(frameMat,nfft));
if odParam.usePowerSpec,
    magSpecMat = magSpecMat .^ 2;
end
diffMagSpecMat = diff([zeros(nfft,1),magSpecMat],1,2);
for i=1:frameNum,
    % calculate spectral flux 
    % (using L1 norm and only considering positive changes)
    diffMagSpec = diffMagSpecMat(:,i);
    SF(i) = sum(diffMagSpec (diffMagSpec > 0) );
end

% step 3. peak picking
% step 3.1 post-processing
% normalize detection function to zero-mean and unit-variance
normalSpecFlux = inputNormalize(SF);

% step 3.2 thresholding 
% (separate event-related and nonevent-related peaks)
% (optional, not finished yet)

% step 3.3 peak-picking 
% (find local maxima)
peak = localMax(normalSpecFlux, odParam.wSize, odParam.mSize, odParam.delta, odParam.alpha);
%decide the time of onsets (divided by 每個frame所佔秒數)
onsetFrameIdx = find(peak==1);
onsetTime = onsetFrameIdx ./ (fs/( frameSize - overlap));
if nargout>2,
    onsetSampleIdx = frame2sampleIndex(onsetFrameIdx, frameSize, overlap);
end

% plotting if necessary
if plotOpt,
    % plot original waveform
    subplot(3,1,1); plot((1:length(y))/fs,y); 
    axis xy; axis tight; ylabel('Amplitude'); title('Waveform');    
  
    % plot spectrogram
    subplot(3,1,2); % spectrogram(y, odParam.frameSize, odParam.frameSize-odParam.overlap, odParam.frameSize, fs, 'yaxis');
    frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
    F = [0:(frameSize/2)]*fs/frameSize;
    imagesc(frameTime,F,20*log10(magSpecMat(1:1+frameSize/2,:)));
    axis xy; ylabel('Frequency(Hz)'); title('FFT Log Magnitude Spectrogram'); colormap(jet); view(0,90);
    
    % plot detection function and detected onsets
    subplot(3,1,3); 
    % plot ground truth, if necessary
    hold on;
    if nargin == 6,
       for i=1:length(groundTruth),
           line(groundTruth(i)*[1 1],[normalSpecFlux(ceil(groundTruth(i)*fs/( frameSize - overlap)+eps)), normalSpecFlux(ceil(groundTruth(i)*fs/( frameSize - overlap)+eps))], 'color', 'k', 'Marker', '*', 'MarkerSize', 8);
       end
    end    
    % plot detected onset
	for i=1:length(onsetTime),
		line(onsetTime(i)*[1 1], [min(normalSpecFlux), max(normalSpecFlux)], 'color', 'r', 'Marker', 'diamond');
    end    
    plot((1:length(SF))/(fs/( frameSize - overlap)), normalSpecFlux);
    axis xy; axis tight; xlabel ('Time(s)'); ylabel('Detection Function'); 
    hold off;
    if nargin < 6,
        title(['Detected ' num2str(length(onsetTime)) ' Onsets (diamond)']);
    else
        title(['Detected ' num2str(length(onsetTime)) ' Onsets (diamond) & Ground Truth ' num2str(length(groundTruth)) ' Onsets(*)']);
    end
    
    % sound for original waveform
	U.y=y; U.fs=fs;
	if max(U.y)>1, U.y=U.y/(2^nbits/2); end  
    % sound for waveform with detected onset information
    % add 5ms noise for each onset time to form a noisy sound
    timeDuration = round(fs*1/1000);    % add 5ms noise
    x = 0:1/timeDuration:2*pi;
    noise = sin(200*x)';
    % hamming windowing
    noise = noise .* hamming(length(noise));
    U.onsetY=y;
    if max(U.onsetY)>1, U.onsetY=U.onsetY/(2^nbits/2); end
    for i=1:length(onsetTime),
        time = ceil(onsetTime(i)*fs+eps);       % +eps ==> avoid zero
        if time>=length(y), continue; end       % the last frame needs not to add noise signals, 20101212
        U.onsetY(time:time+length(noise)-1)=U.onsetY(time:time+length(noise)-1)+noise;
    end
	set(gcf, 'userData', U);
	uicontrol('string', 'Play music', 'callback', 'U=get(gcf, ''userData''); sound(U.y, U.fs);', 'position', [20, 20, 80, 20]);      
    uicontrol('string', 'Play music + detected onsets', 'callback', 'U=get(gcf, ''userData''); sound(U.onsetY, U.fs);', 'position', [120, 20, 180, 20]);    
    
    % sound for waveform with ground truth information, if necessary
    if nargin==6,
        U.truthY=y;
        if max(U.truthY)>1, U.truthY=U.truthY/(2^nbits/2); end
        for i=1:length(groundTruth),
            time = ceil(groundTruth(i)*fs+eps); % +eps ==> avoid zero
            U.truthY(time:time+length(noise)-1) = U.truthY(time:time+length(noise)-1)+noise;
        end
        set(gcf, 'userData', U);
        uicontrol('string', 'Play music + ground-truth onsets', 'callback', 'U=get(gcf, ''userData''); sound(U.truthY, U.fs);', 'position', [320, 20, 200, 20]);
    end
end

% subfunction, find local maxima from detection function
function [peak,cand] = localMax(detFun, w, m, delta, alpha)
frameNum = length(detFun);
% criterion 1, detFun(n) >= detFun(k) for all k (where n-w <= k <= n+w)
detFun1 = [zeros(1,w),detFun,zeros(1,w)];
cand1 = zeros(1,frameNum);
for i=1:frameNum,
    leftFlag = (detFun1(i+w) >= detFun1(i:i+w-1));
    rightFlag = (detFun1(i+w) >= detFun1(i+w+1:i+2*w));
    if  (sum(leftFlag) + sum(rightFlag)) == 2*w,
        cand1(i)=1;
    end
end
% criterion 2, 目前frame的detection value比前m*w+w個values的平均值大超過delta值
detFun2 = [zeros(1,m*w), detFun,zeros(1,w)];
cand2 = zeros(1,frameNum);
for i=1:frameNum,
    tmp = mean(detFun2(i:i+m*w+w));
    if detFun2(i+m*w) >= tmp+delta, cand2(i)=1; end
end
% criterion 3, use a gamma-like funtion to decide a candidate
% PS. see Onset Detection Revisited, Simon Dixon, 2006, for more details
cand3 = zeros(1,frameNum);
g(1)=0;
for i=2:frameNum+1,
    g(i) = max(detFun(i-1),alpha*g(i-1)+(1-alpha)*detFun(i-1));
end
for i=1:frameNum,
    if detFun(i) >= g(i+1), cand3(i)=1; end
end
% decide final candidates
peak = cand1 & cand2 & cand3;
cand=max([cand1;cand2;cand3]);

function selfdemo
waveFile = 'little_star_piano2.wav';
[y,fs,nbits] = wavread(waveFile);
y=y(:,1);
odParam=odParamSet2(fs);
plotOpt=1;

% roundTruth (read cue label)
fid = fopen(waveFile, 'rb');
binWave = fread(fid);
fclose(fid);
cueLabel = binWave2cueLabel(binWave);
groundTruth = cueLabel / fs;

onsetTime=feval(mfilename, y, fs, nbits, odParam, plotOpt, groundTruth);
function spec = NSHS(y, opt, plotOpt)
%NSHS Summary of this function goes here
%   Detailed explanation goes here

if nargin<1, selfdemo; return; end
if nargin<2|isempty(opt)
	opt.fs=16000;
	opt.frameSize=1024;
	opt.overlap=1024-128;
	opt.lowerFreq=80;
	opt.higherFreq=1280;
	opt.weight=0.98;
	opt.fftLen=2^13;
end
if nargin<3, plotOpt=0; end

frameMat=buffer2(y, opt.frameSize, opt.overlap);
frameNum=size(frameMat, 2);
window=hamming(opt.frameSize)*ones(1, frameNum);
frameMat=frameMat.*window;

complexSpec=fft(frameMat, opt.fftLen);		% Complex spectrum
magSpec=abs(complexSpec(1:opt.fftLen/2+1, :));	% Magnitude spectrum

spec=zeros(size(magSpec));
freqVec=(opt.fs/opt.fftLen)*([0:opt.fftLen/2]);				% ©Ò¦³ªº frequency bin
foiIndex=find(opt.lowerFreq<=freqVec & freqVec<=opt.higherFreq);	% Index of FOI (freq. of interest)
index1=min(foiIndex);
index2=max(foiIndex);
weightVec=(opt.weight).^(0:opt.fftLen/2)';

for i=1:frameNum
	for j=index1:index2
		step=j-1;
		toBeAdded=magSpec(j:step:end, i);
		weight=weightVec(1:length(toBeAdded));
		spec(j,i)=sum(toBeAdded.*weight)/sum(weight);
	end
end

if plotOpt==1
	subplot(2,1,1);
	imagesc(db(magSpec)); axis xy;
	title('Original spectrogram');
	xlabel('Frame index');
	ylabel('Frequency bin');
	subplot(2,1,2);
	imagesc(db(spec)); axis xy;
	title('Spectrogram after NSHS');
	xlabel('Frame index');
	ylabel('Frequency bin');
end

% ====== Self demo
function selfdemo
demoFile = 'demo.wav';
y = wavread(demoFile);
spec = nshs(y, [], 1);

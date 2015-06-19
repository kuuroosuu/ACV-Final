function feature=frame2feature(frameMat, fs, opt, plotOpt)
% frame2feature: Frame to feature conversion
%
%	Usage:
%		feature=frame2feature(frameMat, fs, opt, plotOpt)
%
%	Description:
%		feature=frame2feature(frameMat, fs, opt, plotOpt) return the features from a given frame matrix.
%			frameMat: a frame matrix (with each column representing a frame)
%			fs: Sampling rate
%			opt: Options for the feature extraction, whos default value can be obtained via frame2feature('defaultOpt').
%				opt.lowFreq: Lower-bound frequency
%				opt.highFreq: Upper-bound frequency
%				opt.featureType: Type of features to be extracted
%					'specEntropy': Spectrum entropy
%					'specCentroid': Spectrum entropy
%					'specBisectingFreq': Spectrum bisecting frequency
%					'energyInFreqBand': Energy in frequency band
%					'specPeakInFreqBand': Spectral peak energy in frequency band
%					'localMaxInterval': Mean interval between local maxima of the waveform
%					'zeroCrossingRate': Zero crossing rate
%			plotOpt: 1 for plotting
%			feature: Extracted features
%				feature(i).type: Type of feature i (specified by opt.featureType{i})
%				feature(i).value: Value of the extracted feature i
%
%	Example:
%		waveFile='singaporeIsAFinePlace.wav';
%		[y, fs, nbits]=wavread(waveFile);
%		frameSize=640; overlap=480;
%		frameMat=buffer2(y, frameSize, overlap);
%		opt=frame2feature('defaultOpt');
%		feature=frame2feature(frameMat, fs, opt, 1);
%
%	Reference:
%		Jia-lin Shen, Jeih-weih Hung, Lin-shan Lee, "Robust Entropy-based Endpoint Detection for Speech Recognition in Noisy Environments", ICSLP, 1998.

%	Roger Jang, 20120924

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
% ====== Set the default options
if nargin==1 && ischar(frameMat) && strcmpi(frameMat, 'defaultOpt')
	feature.lowFreq=250;
	feature.highFreq=6000;
	feature.featureType={'specEntropy', 'specCentroid', 'specBisectingFreq', 'localMaxInterval'};
	return
end
if nargin<3||isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<4, plotOpt=0; end

if ischar(opt.featureType), opt.featureType={opt.featueType}; end
if size(frameMat,1)==1 || size(frameMat,2)==1, frameMat=frameMat(:); end

[frameSize, frameCount]=size(frameMat);
numFeaType=length(opt.featureType);
% ====== Assume we will extract spectrum-based features...
[magSpec, phaseSpec, freq, powerSpec]=fftOneSide(frameMat, fs);
index=find(opt.lowFreq<freq & freq<opt.highFreq);
freqOrig=freq; freq=freq(index);
magSpecOrig=magSpec; magSpec=magSpec(index, :);
powerSpecOrig=powerSpec; powerSpec=powerSpec(index, :);
totalMagSpec=sum(magSpec);
feature=struct([]);
for j=1:numFeaType
	feature(j).type=opt.featureType{j};
	feature(j).value=zeros(frameCount,1);
	switch(lower(opt.featureType{j}))
		case lower('specEntropy')
			for i=1:frameCount
				p=magSpec(:,i)/totalMagSpec(i);
				feature(j).value(i)=-sum(p.*log2(p));
			end
		case lower('specCentroid')
			feaVec=zeros(frameCount,1);
			for i=1:frameCount
				feature(j).value(i)=dot(magSpec(:,i), freq)/sum(magSpec(:,i));
			end	
		case lower('specBisectingFreq')
			for i=1:frameCount
				cumulated=cumsum(magSpec(:,i));
				index=find(cumulated>totalMagSpec(i)/2);
				feature(j).value(i)=freq(index(1));
			end

		case lower('energyInFreqBand')
			for i=1:frameCount
				feature(j).value(i)=sum(magSpec(:,i).^2);
			end
		case lower('specPeakInFreqBand')
			for i=1:frameCount
				feature(j).value(i)=max(magSpec(:,i));
			end
		case lower('localMaxInterval')
			for i=1:frameCount
				lMaxIndex=find(localMax(frameMat(:,i)));
				feature(j).value(i)=mean(diff(lMaxIndex));
			end
		case lower('zeroCrossingRate')
			for i=1:frameCount
				feature(j).value(i)=sum(frameMat(1:end-1,i).*frameMat(2:end,i)<0);
			end
		otherwise
			error('Unknown feature type %s!\n', opt.featureType);
	end
end

if plotOpt
	plotNum=numFeaType+1;
	subplot(plotNum,1,1); surf(1:frameCount, freqOrig, 10*log10(magSpecOrig), 'edgeColor', 'none'); axis tight; view(0, 90);
	ylabel('Frequency (Hz)'); title('Power spectrum (dB)');
	line([1, frameCount], opt.lowFreq*[1, 1], 'color', 'r');
	line([1, frameCount], opt.highFreq*[1, 1], 'color', 'r');
	for i=2:plotNum
		subplot(plotNum,1,i); plot(feature(i-1).value); title(feature(i-1).type);
		set(gca, 'xlim', [-inf inf]);
	end
	xlabel('Frame index');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);

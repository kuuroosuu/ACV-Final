function ptOpt=ptOptSet(fs, nbits, pfType)
% ptOptSet: Set PT (pitch tracking) parameters
%
%	Usage:
%		ptOpt=ptOptSet(fs, nbits, pfType)
%
%	Description:
%		ptOpt=ptOptSet(fs, nbits, pfType) returns the parameters for pitch tracking.
%			fs: sample rate
%			nbits: bit resolution
%			pfType: 0 for AMDF, 1 for ACF.
%		By default, this function will give a frameSize and overlap that lead to a frameDuration of 32 ms and frameRate of 31.25. You can modify the output as necessary for pitch tracking.
%		If any field names are changed, you also need to update ptByDpOverPfMex.cpp.
%
%	Example:
%		ptOpt=ptOptSet(8000, 8, 1)

if nargin<1, fs=8000; end
if nargin<2, nbits=8; end
if nargin<3, pfType=1; end		% ACF by default, for better performance

%% ====== 固定參數，請勿修改
ptOpt.frameDuration=32;		% Duration (in ms) of a frame
ptOpt.overlapDuration=0;	% Duration (in ms) of overlap
ptOpt.fs=fs;			% 取樣頻率
ptOpt.nbits=nbits;		% 位元解析度
ptOpt.pfType=pfType;		% pitch function type, 0 for AMDF, 1 for ACF
ptOpt.alpha=0.5;		% This is used for pfWeight, used in frame2pitch. Only for ACF and method=1.
				% alpha=0 corresponds exactly to ACF with method=1, which is likely to have double-pitch error
				% alpha=1 corresponds roghtly to ACF with method=2, which is likely to have half-pitch error
				% The best value of alpha is obtained via exhaustive search
ptOpt.usePitchSmooth=1;
ptOpt.resampleRatio=1;		% Resample the waveform. 2 for taking the odd-index samples, etc.
frameRate=31.25;
frameDuration=0.032;		% 32 ms
%% ====== 音框長度及相鄰音框重疊點數
switch(ptOpt.fs)
	case 8000
	%	ptOpt.frameSize=256; ptOpt.overlap=128;	% frame rate = 62.50
		ptOpt.frameSize=256; ptOpt.overlap=0;	% frame rate = 31.25
	case 16000
	%	ptOpt.frameSize=640; ptOpt.overlap=640-160;	% frameSize 是 HTK 的兩倍，frame rate=100, 搭配語音辨識使用
		ptOpt.frameSize=512; ptOpt.overlap=0;	% frame rate = 31.25
	otherwise	% Take care of other fs, assuming frame rate=31.25
		ptOpt.frameSize=round(frameDuration*ptOpt.fs); ptOpt.overlap=ptOpt.frameSize-round(ptOpt.fs/frameRate);
%		if mod(ptOpt.frameSize,2)~=0
%			ptOpt.frameSize=ptOpt.frameSize+1;		% Keep this even for setting pfLen
%		end
end

%% ====== 音高追蹤的一般參數
ptOpt.fs=fs;
ptOpt.nbits=nbits;
ptOpt.pfLen=ptOpt.frameSize/2;						% ACF/AMDF 的平移量
ptOpt.minPitch=30;								% 最低音高
ptOpt.maxPitch=84;								% 最高音高
ptOpt.absVolTh=2*ptOpt.frameSize*(2^(ptOpt.nbits-8));				% 靜音的絕對門檻值
ptOpt.checkMultipleFreq=1;							% 是否進行倍頻轉換
ptOpt.maxAmdfLocalMinCount=20*(ptOpt.frameSize/256)*(ptOpt.fs/8000);		% 若超過此個數，視為氣音（此數需要配合 fs 及 frameSize 調整）
ptOpt.maxPitchDiff=7*(((ptOpt.frameSize-ptOpt.overlap)/ptOpt.fs)/(256/8000));	% 相鄰音高點的最大差值
ptOpt.minAmdfRange=5*(2^nbits/2^8);			% Normalized by nbits
ptOpt.frame2pitchFcn='frame2pitch';			% 由音框計算音高的函數
ptOpt.mainFun='dpOverPfMat';				% Main function for pitch tracking
ptOpt.useSmooth=1;					% Smooth the final output pitch
ptOpt.useRepeat=1;					% Repeat PT twice, using the median of the first for limiting the second PT
ptOpt.useEpd=1;
ptOpt.useEpdBeforePt=0;					% Perform EPD before PT via DP. Valid only when useEpd=1
%% ====== frame normalization
ptOpt.zeroMeanPolyOrder=6;
%% ====== acf/nsdf 的相關參數
ptOpt.clarityTh=0.5;		% k for nsdf
%% ====== dp 的相關參數
ptOpt.dpMethod=1;		% dpMethod=1 ===> 以 indexDiff 來保持平滑度（缺點：對音高是非線性）, dpMethod=2 ===> 以 pitchDiff 來保持平滑度（優點：對音高是線性，缺點：計算量較大）
ptOpt.dpPosDiffOrder=2;		% dpPosDiffOrder=1 ===> abs diff, dpPosDiffOrder=2 ===> squared diff
ptOpt.dpUseLocalOptim=0;	% dpUseLocalOptim=0 ===> use all amdf as ppc, dpUseLocalOptim=1 ===> use only local min of amdf as ppc
ptOpt.ppcNum=10;
ptOpt.pitchRangeMax=24;		% Max range of pitch vector. For singing input
ptOpt.pitchDiffMax=5;		% Max pitch difference between two neighboring frame. For singing input
ptOpt.pfWeight=1;		% Default value, to be modified according to different fs, see below
ptOpt.indexDiffWeight=1;	% Default value, to be modified according to different fs, see below
switch(ptOpt.pfType)
	case 0	% AMDF
		ptOpt.pfWeightedByClarity=1;			% To be verified.
		switch(ptOpt.fs)
			case 8000
				if (ptOpt.frameSize==256) & (ptOpt.overlap==176)
					ptOpt.pfWeight=128;		% 40
					ptOpt.indexDiffWeight=159;	% 540
					ptOpt.dpMethod=1;
					ptOpt.dpPosDiffOrder=1;
					ptOpt.dpUseLocalOptim=0;
				end
				if (ptOpt.frameSize==256) & (ptOpt.overlap==0)	% Obtained from goDpParamOptim.m, frameRate=31.25
					ptOpt.pfWeight=1;
					ptOpt.indexDiffWeight=15;
					ptOpt.dpMethod=1;
					ptOpt.dpPosDiffOrder=1;
					ptOpt.dpUseLocalOptim=0;
				end
				if (ptOpt.frameSize==256) & (ptOpt.overlap==128)	% Obtained from goDpParamOptim.m, frameRate=62.5
					ptOpt.pfWeight=100;
					ptOpt.indexDiffWeight=36;
					ptOpt.dpMethod=1;
					ptOpt.dpPosDiffOrder=2;
					ptOpt.dpUseLocalOptim=0;
				end
			case 16000
				if (ptOpt.frameSize==512) & (ptOpt.overlap==0)
					% ====== Parameters for pitchDiffSquare
					ptOpt.pfWeight=32;
					ptOpt.indexDiffWeight=321;
				end
				if ((ptOpt.frameSize==640) & (ptOpt.overlap==480))
					ptOpt.pfWeight=12800;
					ptOpt.indexDiffWeight=700;	% original: 700
					ptOpt.dpMethod=2;
					ptOpt.dpPosDiffOrder=2;
					ptOpt.dpUseLocalOptim=1;
				end
			otherwise
			%	fprintf('Using default parameters (which might not be complete) for fs=%d\n', ptOpt.fs);
		end
	case 1	% ACF
		ptOpt.pfWeightedByClarity=1;			% based on a limited test
		switch(ptOpt.fs)
			case 8000
				if (ptOpt.frameSize==256) & (ptOpt.overlap==0)
					ptOpt.pfWeight=4;
					ptOpt.indexDiffWeight=388;
					ptOpt.dpMethod=2;
					ptOpt.dpPosDiffOrder=1;
					ptOpt.dpUseLocalOptim=1;
					
					ptOpt.pfWeight=250;
					ptOpt.indexDiffWeight=833;
					ptOpt.dpMethod=2;
					ptOpt.dpPosDiffOrder=2;
					ptOpt.dpUseLocalOptim=0;
				end
				if (ptOpt.frameSize==256) & (ptOpt.overlap==128)	% Obtained from goDpParamOptim.m, frameRate=62.5
					ptOpt.pfWeight=1;
					ptOpt.indexDiffWeight=80;
					ptOpt.dpMethod=2;
					ptOpt.dpPosDiffOrder=1;
					ptOpt.dpUseLocalOptim=1;
				end
				if (ptOpt.frameSize==256) & (ptOpt.overlap==176)
					ptOpt.pfWeight=1;
					ptOpt.indexDiffWeight=345;
					ptOpt.dpMethod=2;
					ptOpt.dpPosDiffOrder=1;
					ptOpt.dpUseLocalOptim=0;
				end
			case 16000
				if (ptOpt.frameSize==640) & (ptOpt.overlap==480)
				%	ptOpt.pfWeight=1;
				%	ptOpt.indexDiffWeight=4;
				%	ptOpt.dpMethod=2;
				%	ptOpt.dpPosDiffOrder=2;
				%	ptOpt.dpUseLocalOptim=1;
					ptOpt.pfWeight=1;
					ptOpt.indexDiffWeight=6;
					ptOpt.dpMethod=1;
					ptOpt.dpPosDiffOrder=2;
					ptOpt.dpUseLocalOptim=0;
				end
			otherwise
			%	fprintf('Using default parameters (which might not be complete) for fs=%d\n', ptOpt.fs);
		end
	case 2	% NSDF
		% Nothing here
	case 3	% SHS
		% Nothing here	
	otherwise
		error(sprintf('Unsupported pfType=%d!', ptOpt.pfType));
end
%% ====== Parameters for removing unreliable pitch (due to silence or unvoiced sounds)
ptOpt.useHmm4suvDetection=0;
ptOpt.useVolThreshold=1;		% volRatio is determine by epdByVol.
ptOpt.useClarityThreshold=1;
ptOpt.clarityRatio=0.5;			% Effective only when useClarityThreshold=1
ptOpt.localMaxCountThresholding=0;
ptOpt.localMaxCountRatio=0.5;		% Effective only when localMaxCountThresholding=1
%% ====== Parameters for pre processing
ptOpt.useWaveEnhancement=0;		% Apply wave enhancement
ptOpt.useLowPassFilter=0;		% Apply low-pass filtering on the input waveform
ptOpt.useHighPassFilter=0;		% Apply high-pass filtering on the input waveform
%% ====== Parameters for post processing
ptOpt.medianFilterOrder=5;		% Use median filter to smooth pitch (Not necessary for pitchTrackingForcedSmooth.m)
ptOpt.useParabolicFit=1;		% Apply parabola for finding the min/max of the pitch function
ptOpt.cutLeadingTrailingZero=0;
%% ====== Parameter for forcing smooth by changing indexDiffWeight (used in pitchTrackingForcedSmooth.m)
ptOpt.pitchDiffTh=8;
ptOpt.roundNum=12;
ptOpt.searchMethod='binarySearch';	% 'binarySearch' or 'linearSearch'
ptOpt.leftBound=0;
ptOpt.rightBound=5000;	% for binary search only
ptOpt.boundRangeTh=10;
%% ====== Parameters for PT performance evaluation
ptOpt.pitchTol=1;			% Tolerance for computing pitch accuracy
ptOpt.useSil4perfEval=0;
%% ====== Derived quantities
ptOpt.frameRate=ptOpt.fs/(ptOpt.frameSize-ptOpt.overlap);
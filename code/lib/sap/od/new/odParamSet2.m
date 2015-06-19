function odParam=odParamSet2(fs)
% odParamSet: Set parameters for onset detection
% 	Usage: od=mfccParamSet(fs)

%	Jia-Min Ren, 20091101, 20091124

if nargin<1; fs=16000; end

odParam.frameSize=round(fs/(44100/2048));	% Frame size (fs=44100 ===> frameSize=2048, 46ms)
odParam.overlap=round(fs/(44100/1637));     % overlap, 36ms
odParam.nfft = odParam.frameSize;           % number of fft bin
odParam.usePowerSpec=0;                        % use power spectrum or not, 0:no, 1:yes

% The followings are mainly for onsetDetection.m
odParam.downSample=0;            % 0: do not downsample; 
odParam.nmel = 40;                  % mel channels

% The followings are mainly for peak-picking stage in odBySpecFlux.m and
% odByPhaseDev.m
odParam.wSize = 3;          % # of window used to find a local maximum
odParam.mSize = 3;          % calculate the mean over a larger range before the peak
odParam.delta = 1;          % above the local mean which an onset must reach
odParam.alpha = 0.7;        % the weight of considering current detection value or previous one

% The followings are mainly for odByPhaseDev.m
odParam.phaseDev = 0;       % 0 for phase deviation
                            % 1 for weighted phase deviation
                            % 2 for normalized weighted phase deviation
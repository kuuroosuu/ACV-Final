function odPrm=odPrmSet(fs)
% odPrmSet: Set parameters for onset detection
% 	Usage: od=mfccPrmSet(fs)

%	Jia-Min Ren, 20091101

if nargin<1; fs=16000; end

odPrm.frameSize=round(fs/(44100/2048));	% Frame size (fs=44100 ===> frameSize=2048, 46ms)
odPrm.overlap=round(fs/(44100/1637));     % overlap, 36ms
odPrm.nfft = odPrm.frameSize;           % number of fft bin

% The followings are mainly for onsetDetection.m
odPrm.downSample=0;            % 0: do not downsample; 
odPrm.nmel = 40;                  % mel channels

% The followings are mainly for peak-picking stage in odBySpecFlux.m and
% odByPhaseDev.m
odPrm.wSize = 3;          % # of window used to find a local maximum
odPrm.mSize = 3;          % calculate the mean over a larger range before the peak
odPrm.delta = 2.25;          % above the local mean which an onset must reach
odPrm.alpha = 0.5;        % the weight of considering current detection value or previous one

% The followings are mainly for odByPhaseDev.m
odPrm.phaseDev = 2;       % 0 for phase deviation
                            % 1 for weighted phase deviation
                            % 2 for normalized weighted phase deviation
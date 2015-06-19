function btOpt=btOptSet

addpath /users/jang/matlab/toolbox/utility
addpath /users/jang/matlab/toolbox/machineLearning
addpath /users/jang/matlab/toolbox/sap

btOpt.type='time-varying';	% 'constant' or 'time-varying' tempo
btOpt.wingRatio=0.1;		% One-side tolerance for beat position identification via forward/forward search (only for constant tempo)
btOpt.bpmMax=250;
btOpt.bpmMin=50;		% 70
btOpt.trialNum=8;		% No. of trials for beat positions
btOpt.acfMethod=2;		% Method for computing NC's ACF
btOpt.useDoubleBeatConvert=0;	% Double beat conversion: Convert to double beat if necessary
btOpt.useTripleBeatConvert=0;	% Triple beat conversion: Convert to double beat if necessary
btOpt.peakHeightTol=0.11;	% If two peaks are different by 10% with multiple BPM, small BPM is selected.

btOpt.waveDir='D:\dataSet\beatTracking\mirex-beatTracking\mirex06train';
btOpt.waveDir='D:\dataSet\beatTracking\BeatlesWav';
btOpt.waveDir='D:\users\jang\app\beatTracking\BeatlesShortClip';
btOpt.waveDir='D:\users\jang\app\beatTracking\BeatlesShortClip4training';
btOpt.DP_penalty=0.01;	% Transition penalty used in the DP for finding the tempo curve (default:0.01)
btOpt.alpha=1.7;	% Controls how tightly the estimated tempo is enforced within the search range during beat tracking. (default:1.7)  
btOpt.useAcf=0;		% Determines that the Fourier tempogram will combine with the auto-correlation tempogram or not. (0:No 1:Yes)

[parentDir, mainName]=fileparts(btOpt.waveDir);
btOpt.outDir=sprintf('%s_output', mainName);
btOpt.waveDataFile=[btOpt.outDir, '/waveData'];

btOpt.tolerance=0.07;
btOpt.FFT_mode=1;	% 0:Full FFT 1:Partial FFT
btOpt.oscOpt=wave2osc('defaultOpt');				% Options for onset strength curve

btOpt.weight_window_len = 8/2;	%!!!
%sgsrate = sro/shop; % sampling rate of onset strength waveform (onsetstrwaveform)
btOpt.sgsrate = btOpt.oscOpt.fsTarget/btOpt.oscOpt.frameStep/2;		%!!!		% Sample rate of the novelty curve

btOpt.ncFrameSize=btOpt.weight_window_len*btOpt.sgsrate;	% 512
btOpt.ncFrameStep=8;						% overlap=504

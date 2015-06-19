function [tempogramObj, tempogramLocalMax]=osc2tempogram(osc, btOpt, plotOpt)
% osc2tempogram: Compute tempogram from a given onset strength curve (OSC) 
%
%	Usage:
%		[tempogramObj, tempogramLocalMax]=osc2tempogram(osc, btOpt, plotOpt)
%
%	Example:
%	%	waveFile='\dataSet\beatTracking\mirex06train\train1.wav';
%		waveFile='song01s5.wav';
%		btOpt=btOptSet;
%		[ncObj, ncTime]=wave2osc(waveFile, btOpt.oscOpt);
%		tempogram=osc2tempogram(ncObj, btOpt, 1);

%	Roger Jang, 20120330

if nargin<1, selfdemo; return; end
if nargin<2, btOpt=btOptSet; end
if nargin<3, plotOpt=0; end

weight_window_len=btOpt.weight_window_len;
ncFs=btOpt.sgsrate;
ncFrameSize=btOpt.ncFrameSize;	% 512
ncOverlap=btOpt.ncFrameSize-btOpt.ncFrameStep;
% ======================== Fourier tempogram estimation ===========================
zpad = round(weight_window_len/2 * ncFs);
Zpaded_onset_str = [zeros(1,zpad),osc.signal,zeros(1,zpad)];
[tempogram0,ff,tempogramObj.time_index] = mySpecgram(Zpaded_onset_str,ncFrameSize,ncFs,ncFrameSize,ncOverlap, btOpt.FFT_mode, true);
    
switch btOpt.FFT_mode	% 0:Full FFT 1:Partial FFT
    case 1		% Partial FFT
        freq_step = ncFs/(32*ncFrameSize);
        tempo_start = 0.5/freq_step;
        tempo = ff(1+tempo_start:end) * 60;
        rTempogram = tempogram0(1+tempo_start:end,:);
    case 0		% Full FFT
        freq_step = ncFs/ncFrameSize;
        tempo_start = 0.5/tempogramObj;
        tempo_end = 8.2/freq_step;
        tempo = ff(1+tempo_start:tempo_end) * 60;
        rTempogram = tempogram0(1+tempo_start:tempo_end,:);
    otherwise
        error('FFT Mode setting error (0:Full FFT 1:Partial FFT)');
end
tempogramObj.time_index=tempogramObj.time_index/2;	%!!!
tempogramObj.tempo_index = tempo;
tempogramObj.Tempogram = abs(rTempogram);

% ================ Combine Auto-correlation tempogram (optional) ==================
if (btOpt.useAcf)
    [Acorr_tempocurve Acorr_tempoindex] = Autocorr_tempocurve(osc.signal,[30 500],fs);
    Autocorr_TempoCurve_FitFourier = interp1(Acorr_tempoindex,Acorr_tempocurve,tempogramObj.tempo_index);
    Autocorr_TempoCurve_FitFourier = Autocorr_TempoCurve_FitFourier';
    tempogramObj.Tempogram = tempogramObj.Tempogram .* ...
        repmat(Autocorr_TempoCurve_FitFourier,1,length(tempogramObj.time_index));
end

% =========================== Tempo curve estimation ==============================
tempogramLocalMax = localMax(tempogramObj.Tempogram);	% Peak Picking

if plotOpt
	figure;
	subplot(211);
	imagesc(tempogramObj.time_index, tempo, tempogramObj.Tempogram);
	colormap(gray); axis xy
	xlabel('Time (sec)'); ylabel('Beats Per Minute (BPM)'); title('Tempogram');
	subplot(212); 
	[nr,nc] = size(tempogramLocalMax);
	Matrix_tempo=repmat(tempo, nc, 1);
	Peak_tempo_temp=Matrix_tempo.*tempogramLocalMax';
	Peak_tempo=Peak_tempo_temp';
	for i = 1:nr
		line(tempogramObj.time_index, Peak_tempo(i,:), 'marker', '.', 'linestyle', 'none', 'markersize', 1);
	end
	axis([-inf, inf, 0, 500]); axis xy
	xlabel('Time(sec)'); ylabel('Tempo (BPM)'); title('Local maxima of tempogram');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);

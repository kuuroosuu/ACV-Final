function beatTime=tempoCurve2beat(tempoCurve, oscObj, tempogramObj, btOpt, plotOpt)
% tempoCurve2beat:
%
%	Usage:
%		beatTime=tempoCurve2beat(tempoCurve, oscObj, tempogramObj, btOpt, plotOpt)
%
%	Example:
%		waveFile='D:\dataSet\beatTracking\mirex06train\train1.wav';
%		btOpt=btOptSet;
%		oscObj=wave2osc(waveFile, btOpt.oscOpt);
%		tempogram=osc2tempogram(oscObj, btOpt);
%		tempoCurve=tempogram2tempoCurve(tempogram, btOpt);
%		gtBeatFile=strrep(waveFile, 'wav', 'txt');
%		btOpt.gtBeat=gtBeatRead(gtBeatFile);
%		beatTime=tempoCurve2beat(tempoCurve, oscObj, tempogram, btOpt, 1);

if nargin<1, selfdemo; return; end
if nargin<5, btOpt=btOptSet; end
if nargin<6, plotOpt=0; end

% ================================ Beat tracking ==================================
% Dynamic Programming on Beat Tracking
Peak_onset = localMax(oscObj.signal);
Peak_onset_indx = find(Peak_onset == true);
Peak_str = oscObj.signal(Peak_onset_indx);
[dummy,sorted_indx] = sort(Peak_str,'descend');
DP_begin_indx = Peak_onset_indx(sorted_indx(2));
clear Peak_onset Peak_onset_indx Peak_str Peak_val max_indx;

% !!!
Backtrack_beats_backward_indx = beatsearchDP(oscObj.signal, DP_begin_indx, tempoCurve, btOpt.sgsrate*2,tempogramObj.time_index,2,[0.4 1.9],btOpt.alpha,1);
Backtrack_beats_forward_indx =  beatsearchDP(oscObj.signal, DP_begin_indx, tempoCurve, btOpt.sgsrate*2,tempogramObj.time_index,2,[0.4 1.9],btOpt.alpha,2);

if isempty(Backtrack_beats_backward_indx)
    Backtrack_beats_indx = Backtrack_beats_forward_indx;
else if isempty(Backtrack_beats_forward_indx)
        Backtrack_beats_indx = fliplr(Backtrack_beats_backward_indx);
    else
        Backtrack_beats_indx = [fliplr(Backtrack_beats_backward_indx),Backtrack_beats_forward_indx(2:end)];
    end
end
beatTime = oscObj.time(Backtrack_beats_indx);
% ====== Re-estimate the tempo
ibi=diff(beatTime);		% Inter-beat interval
reEstimatedTempo=60./ibi;			% BPM
reEstimatedTempoTime=(beatTime(1:end-1)+beatTime(2:end))/2;

if plotOpt
	figure;
	subplot(211);
	plot(tempogramObj.time_index, tempoCurve, '.-');
	line(reEstimatedTempoTime, reEstimatedTempo, 'marker', 'o', 'color', 'r');
	xlabel('Time (sec)'); ylabel('Tempo (BPM)'); title('Tempo curve');
	legend('Tempo curve', 'Re-estimated tempo curve');
	subplot(212);
	plot(oscObj.time, oscObj.signal); 
	line(beatTime, oscObj.signal(Backtrack_beats_indx), 'marker', 'o', 'color', 'r', 'linestyle', 'none');
	set(gca, 'xlim', [-inf inf]);
	xlabel('Time (sec)'); ylabel('Amplitude'); title('Computed beat positions');
	legend('Novelty curve', 'Identified beat positions');
	if isfield(btOpt, 'gtBeat')
		axisLimit=axis;
	%	for i=1:length(btOpt.gtBeat)
		for i=1:1
			for j=1:length(btOpt.gtBeat{i})
				line(btOpt.gtBeat{i}(j)*[1 1], axisLimit(3:4), 'color', getColorLight(i));
			end
		end
		legend('Novelty curve', 'Identified beats', 'Desired beats');
	end
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);

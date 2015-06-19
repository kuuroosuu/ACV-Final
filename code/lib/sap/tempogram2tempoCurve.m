function tempoCurve=tempogram2tempoCurve(tempogramObj, btOpt, plotOpt)

if nargin<1, selfdemo; return; end
if nargin<2, btOpt=btOptSet; end
if nargin<3, plotOpt=0; end

% Dynamic Programming
tempogramLocalMax = localmax(tempogramObj.Tempogram,0);	% Peak Picking
Peak_data = sparse(tempogramLocalMax.*tempogramObj.Tempogram);
Backtrack = DP(Peak_data, btOpt.DP_penalty, 2);
tempoCurve = tempogramObj.tempo_index(Backtrack);

if plotOpt
	figure;
	plot(tempogramObj.time_index, tempoCurve, '.-');
	BacktrackSmooth = DP(Peak_data, 0.03, 2);
	tempoCurve2 = tempogramObj.tempo_index(BacktrackSmooth);
	line(tempogramObj.time_index, tempoCurve2, 'marker', 'o', 'color', 'r', 'linestyle', 'none');
	xlabel('Time (sec)'); ylabel('Tempo (BPM)'); title('Tempo curves');
	legend('\Theta = 0.01','\Theta = 0.03');
end

% ====== Self demo
function selfdemo
waveFile='D:\dataSet\beatTracking\mirex06train\train1.wav';
btOpt=btOptSet;
[nc, ncTime]=wave2osc(waveFile, btOpt.oscOpt);
tempogram=osc2tempogram(nc, btOpt);
tempoCurve=tempogram2tempoCurve(tempogram, btOpt, 1);

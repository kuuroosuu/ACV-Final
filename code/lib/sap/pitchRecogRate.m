function [recogRate, absAveError, pvData]=pitchRecogRate(pvData, ptOpt)
% pitchRecogRate: PT on a set of wav/pv files and compute the performance
%	Usage: [recogRate, absAveError, pvData]=pitchRecogRate(pvFiles, ptOpt)

%	Roger Jang, 20060214, 20070216

for i=1:length(pvData)
	waveFile=[pvData(i).path(1:end-2), 'wav'];
	wObj=waveFile2obj(waveFile);
%	wObj.signal=wObj.signal(1:2:end); wObj.fs=8000;		% Convert to 8k16b
	t0=clock;
	ptOpt.tPitch=pvData(i).tPitch;	% This is for plotting the target pitch in the next line
	pvData(i).cPitch=feval(ptOpt.method, wObj, ptOpt, 1);		% Pitch tracking
%	fprintf('\tPress key to continue...'); pause; fprintf('\n');
	pvData(i).time=etime(clock, t0);
	[pvData(i).recogRate, pvData(i).absAveError, pvData(i).correctCount, pvData(i).pitchCount]=pvSimilarity(pvData(i).cPitch, pvData(i).tPitch, ptOpt.pitchTol);
%	fprintf('\t%d/%d: (%g sec, %.2f%%) %s\n', i, length(pvData), pvData(i).time, pvData(i).recogRate*100, pvData(i).path);
end

absAveError=sum([pvData.absAveError].*[pvData.pitchCount])/sum([pvData.pitchCount]);
recogRate=sum([pvData.correctCount])/sum([pvData.pitchCount]);
fprintf('\t\tTotal time = %.2f ¬í\n', sum([pvData.time]));
fprintf('\t\tAverage time = %.2f ¬í\n', sum([pvData.time])/length(pvData));
fprintf('\t\tAverage recog. rate = %.2f%%\n', 100*recogRate);
addMyPath;

waveFile='yankeeDoodle_8k8b.wav';
waveFile='ivegotrhythm_16k16b.wav';
waveFile='..\山月照彈琴.wav';
waveFile='..\但使龍城飛將在.wav';
fprintf('waveFile=%s\n', waveFile);

[y, fs, nbits]=wavReadInt(waveFile);
PP=ptOptSet(fs, nbits);
framedY=buffer(y, PP.frameSize, PP.overlap);
frameNum=size(framedY,2);
for i=1:frameNum
%	fprintf('frameIndex=%d/%d\n', i, frameNum);
	frame=framedY(:, i);
	frame2=frameFlip(frame);
%	frame2=localAverage(frame2);
%	amdf=round(frame2amdf(frame2, PP.maxShift, 2, 0));	% Mimic MCU
	amdf=frame2smdf(frame2, PP.maxShift, 3);
	index1=amdf2ppc(   amdf, fs, nbits, PP);
	index2=amdf2ppcMex(amdf, fs, nbits);
	value1=amdf(index1);
	value2=amdf(index2);
%	plot(1:PP.maxShift, amdf, '.-', index1, value1, '^r', index2, value2, 'V');
%	legend('amdf', 'ppc by amdf2ppc', 'ppc by amdf2ppcMex');
	if ~isequal(value1, value2)
		fprintf('Errors!\n');
		fprintf('index1=%s\n', mat2str(index1));
		fprintf('index2=%s\n', mat2str(index2));
		fprintf('value1=%s\n', mat2str(value1));
		fprintf('value2=%s\n', mat2str(value2));
		fprintf('frameIndex=%d\n', i);
		return
	end
end
fprintf('Test done on %d frames, no problem is found!\n', frameNum);
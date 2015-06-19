function pv2midiFile(pitchVec, pitchRate, midiFile)
% pv2midiFile: Pitch vector to midi file conversion
%	Usage: pv2midiFile(pitchVec, pitchRate, midiFile)
%		pitchVec: pitch vector
%		pitchRate: number of pitch points per second (which is equal to fs/(frameSize-overlap))
%		midiFile: output midi file
    
note = pv2note(pitchVec);	% Note segmentation. The format of noteVec is [pitch, duration, pitch, duration...]
noteVec = zeros(1,2*length(note));
for i = 1 : length(note)
    noteVec((i-1)*2+1) = note(i).pitch;
    noteVec(i*2) = note(i).duration;
end

if noteVec(1)==0
	noteVec = noteVec(3:end);	% Get rid of leading rest
end
if noteVec(end-1) == 0
	noteVec = noteVec(1:end-2);	% Get rid of trailing rest
end
    
for i=1:2:length(noteVec)
	k=(i+1)/2;				% Note index
	pVector(k)=noteVec(i);			% pitch
	durVector(k)=noteVec(i+1)*pitchRate;	% duration
	if k==1
		onVector(k)=0;			% onset begin at 0
	else
		onVector(k)=onVector(k-1)+durVector(k-1);
	end
end

%pVector=noteVec(1:2:end);
%durVector=noteVec(2:2:end);
%onVector=[0, cumsum(durVector)];
    
beatUnit = 2*mean(durVector);	% 取音長平均值的兩倍當作一拍
nmat=[];			% nmat=[beatOn beatDur channel pitch velocity beatOnSec beatDurSec ]
 
for i=1:length(pVector)
	beatOn = onVector(i)/beatUnit;		% OnSet 的 beat
	beatDur = durVector(i)/beatUnit;	% 音長的 beats
	channel = 1;				% 第一軌
	vel = 110;				% 音量固定110
	if pVector(i)~= 0 
		nmat=[nmat; beatOn beatDur channel round(pVector(i)) vel onVector(i) durVector(i)];
	end
end

lastOn = nmat(end,1)+nmat(end,2);	% 結尾的拍數
lastSec = nmat(end,6)+nmat(end,7);	% 結尾的時間
tempo = floor(lastOn*60/lastSec);	% 計算節拍速度
writemidi(nmat, midiFile, 120, tempo);	% 寫入檔案，使用miditoolbox
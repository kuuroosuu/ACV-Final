clear all; close all
addpath ../spce4singing

waveFile='../spce4singing/waveData/cweb/frida_faraway/this old man.wav';
waveFile='../spce4singing/waveData/cweb/roger_faraway/the more we get together.wav';
waveFile='../spce4singing/waveData/cweb/beball/the more we get together.wav';
waveFile='../spce4singing/waveData/cweb/frida/the more we get together_ta.wav';
waveFile='../spce4singing/waveData/cweb/roger/the more we get together.wav';
%waveFile='../spce4speech/waveData/frida/Does it make sense to let children play with matches_70.wav';
[PP, CP]=setPrm;
PP.duration=15;
plotOpt=1;
[pitch0, volume, PP] = wave2pitchVolume(waveFile, PP.duration, plotOpt, PP);

[y, fs, nbits]=waveFileRead(waveFile);
framedY=buffer2(y, PP.frameSize, PP.overlap);
frameNum=size(framedY,2);
amdfMat=[];
probMap=[];
ppcNum=5;

% ====== 計算每一個音框的 PPC (Pitch Period Candidate)
for i=1:frameNum
	frame=framedY(:,i);
	volume(i)=frame2volume(frame);
	frame2=frameFlip(frame);
	frame3=localAverage(frame2);
	amdf=frame2amdf(frame3);
	ppcIndex=amdf2ppc(amdf, ppcNum);
	for j=1:length(ppcIndex)
		A(i).ppc(j).index=ppcIndex(j);
		A(i).ppc(j).amdfValue=amdf(ppcIndex(j));
		A(i).ppc(j).pitch=10*freq2pitch(fs/(A(i).ppc(j).index-1));
	end
	amdfMat=[amdfMat, amdf];
end
%volumeThreshold=getVolumeThreshold(volume);
%lowVolIndex=find(volume<volumeThreshold);
%A(lowVolIndex)=[];
%pitch0(lowVolIndex)=[];
%amdfMat(:, lowVolIndex)=[];
%probMap(:, lowVolIndex)=[];

% ====== 開始進行 dynamic programming
w=100;
fprintf('Computing DP table...\n');
% 計算第一列
i=1;
for j=1:length(A(i).ppc)
	A(i).ppc(j).totalDist=0;
	A(i).ppc(j).from=[];
end
% 計算其它各列
for i=2:frameNum
	for j=1:length(A(i).ppc)
		minDist=inf;
		for k=1:length(A(i-1).ppc)
			thisDist=A(i-1).ppc(k).totalDist+w*(A(i-1).ppc(k).index-A(i).ppc(j).index)+A(i).ppc(j).amdfValue;
			if thisDist<minDist
				minDist=thisDist;
				A(i).ppc(j).from=[i-1, k]';
			end
		end
		A(i).ppc(j).totalDist=minDist;
	end
end

% 找出最後一列的最小值
[minDist, minIndex]=min([A(end).ppc.totalDist]);
% Back tracking 以找出最佳路徑
i=frameNum;
j=minIndex;
path=[i; j];
while ~isempty(A(i).ppc(j).from)
	from=A(i).ppc(j).from;
	path=[from, path];
	i=from(1);
	j=from(2);
end

% 查出最佳路徑的 pitch
pitch=[];
for i=1:size(path,2)
	p=path(1,i);
	q=path(2,i);
	pitch=[pitch, A(p).ppc(q).pitch];
end
figure; plot(1:length(pitch0), pitch0, '.-', 1:length(pitch), pitch, 'o-');
legend('wave2pitchVolume', 'DP');

return

% 找出最佳路徑的座標
x=[];
y=[];
for i=1:size(path,2)
	p=path(1,i);
	q=A(i).ppc(path(2,i)).index;
	x=[x, p];
	y=[y, q];
end

% 畫出 amdfMat
figure; imagesc(amdfMat); colorbar; axis xy; title('amdfMat');
line(x, y, 'color', 'w', 'marker', 'o');
% 畫出所有的 local minima
for i=1:length(A)
	for j=1:length(A(i).ppc)
		line(i, A(i).ppc(j).index, 'color', 'k', 'marker', '.', 'linestyle', 'none'); 
	end
end
% 畫出 probMap w.r.t. index
figure; imagesc(probMap); colorbar; axis xy; title('probMap w.r.t. index');
line(x, y, 'color', 'w', 'marker', 'o');
% 畫出所有的 local minima
for i=1:length(A)
	for j=1:length(A(i).ppc)
		line(i, A(i).ppc(j).index, 'color', 'k', 'marker', '.', 'linestyle', 'none'); 
	end
end
% 畫出 probMap w.r.t. pitch
pitchIndex=1:size(probMap,1);
frameIndex=1:size(probMap,2);
pitchIndex=freq2semitone(8000./(pitchIndex-1));
figure;
surf(pitchIndex, frameIndex, probMap'); view([-90, 90]); shading flat; title('probMap w.r.t. pitch');
xlabel('Frame index'); ylabel('Pitch');
return
line(x, y, 'color', 'w', 'marker', 'o');
% 畫出所有的 local minima
for i=1:length(A)
	for j=1:length(A(i).ppc)
		line(i, A(i).ppc(j).index, 'color', 'k', 'marker', '.', 'linestyle', 'none'); 
	end
end
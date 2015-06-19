function [pitch, A]=A2pitchByDp(A, fs, ptOpt)
% vowel2pitchByDp: 利用 DP 的方法來找出最佳的音高曲線
%	Usage: [pitch, A]=A2pitchByDp(A, fs, ptOpt)

%	Roger Jang, 20050228

if nargin<2, fs=8000; end
if nargin<3, ptOpt=ptOptSet(fs, nbits); end

% Generate the table from ppcIndex to pitch (which is integer)
% This table is used in DP
index2pitch=round(10*freq2pitch(ptOpt.fs./(1:ptOpt.frameSize/2-1)));		% 乘以 10 以後的整數音高值！
index2pitch=[0, index2pitch];

frameNum=length(A);

% ====== 根據 A 開始進行 dynamic programming
% ====== 第一列
ppcNum=length(A(1).ppc);
for j=1:ppcNum
	A(1).ppc(j).totalDist = ptOpt.amdfWeight*A(1).ppc(j).amdfValue;
	A(1).ppc(j).from = -1;
end
% ====== 其他列	
for i=2:frameNum
%	fprintf('i=%d/%d\n', i, frameNum);
	for j=1:length(A(i).ppc)
		minDist = inf;
		for k=1:length(A(i-1).ppc)
			if ptOpt.dpMethod==1
				indexDiff=abs(A(i-1).ppc(k).index-A(i).ppc(j).index);
				if ptOpt.dpPosDiffOrder==1
					thisDist = A(i-1).ppc(k).totalDist+ptOpt.indexDiffWeight*(indexDiff);
				else
					thisDist = A(i-1).ppc(k).totalDist+ptOpt.indexDiffWeight*(indexDiff*indexDiff);
				end
			else
				pitchDiff=abs(index2pitch(A(i-1).ppc(k).index)-index2pitch(A(i).ppc(j).index));
				if ptOpt.dpPosDiffOrder==1
					thisDist = A(i-1).ppc(k).totalDist+ptOpt.indexDiffWeight*(pitchDiff);
				else
					thisDist = A(i-1).ppc(k).totalDist+ptOpt.indexDiffWeight*(pitchDiff*pitchDiff);
				end
			end
			if thisDist < minDist
				minDist = thisDist;
				A(i).ppc(j).from = k;
			end
		end
		A(i).ppc(j).totalDist = minDist + ptOpt.amdfWeight*A(i).ppc(j).amdfValue;

%		if j==2
%			keyboard
%		end

	end
end
	
% ====== 找出最後一列的最小值
minDist = inf;
i=frameNum;
for j = 1:length(A(i).ppc)
	if A(i).ppc(j).totalDist<minDist
		minDist = A(i).ppc(j).totalDist;
		minIndex = j;
	end
end
	
% ====== 查出最佳路徑的 pitch period
for i=frameNum:-1:1
	A(i).pp = A(i).ppc(minIndex).index;
	minIndex = A(i).ppc(minIndex).from;
end

% ====== 轉成 semitone
for i = 1:frameNum
	pitch(i) = 69+12*log2(fs/(A(i).pp-1)/440);
end
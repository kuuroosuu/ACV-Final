function minDist=dtw4qbsh(pitch, song, qbshPrm)
%dtw4qbsh: DTW for QBSH
%	Usage: minDist=dtw4qbsh(pitch, song, qbshPrm)
%		pitch: input pitch vector
%		song: pitch vector of a song in the song database

%	Roger Jang, 20090922

pitch=pitch(:)';
song=song(:)';
%% ====== Create options for DTW
dtwOpt=dtw('defaultOpt');
dtwOpt.beginCorner=qbshPrm.beginCorner;
dtwOpt.endCorner=qbshPrm.endCorner;
switch qbshPrm.method
	case 'dtw1'
		dtwOpt.type=1;
	case 'dtw2'
		dtwOpt.type=2;
	otherwise
		error('Unknown dtw function!\n');
end
%% ====== Invoke dtwCore4qbsh()
switch qbshPrm.matchType
	case 'wave2midi'
		minDist=dtwCore4qbsh(pitch, song, dtwOpt, qbshPrm.dtwCount);
	case 'wave2wave'
		minDist1=dtwCore4qbsh(pitch, song, dtwOpt, qbshPrm.dtwCount);
		minDist2=dtwCore4qbsh(song, pitch, dtwOpt, qbshPrm.dtwCount);
		minDist=min(minDist1, minDist2);
	otherwise
		error('Unknwon matchType!');
end

% === DTW core function for QBSH with key transposition via binary-like search
function minDist=dtwCore4qbsh(pitch, song, dtwOpt, dtwCount)
searchRange=1;		% search range = plus/minus 2
dist1=dtw(pitch, song, dtwOpt);
for j=1:(dtwCount-1)/2
	dist2 = dtw(pitch-searchRange, song, dtwOpt);
	dist3 = dtw(pitch+searchRange, song, dtwOpt);
	dist=[dist1, dist2, dist3];
	shift=[0, -searchRange, searchRange];
	[dist1, minIndex]=min(dist);
	pitch=pitch+shift(minIndex);
	searchRange=searchRange/2;
end
minDist=dist1;
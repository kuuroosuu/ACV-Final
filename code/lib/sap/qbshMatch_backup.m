function [distVec, scoreVec] = qbshMatch(songData, pitch, mrParam)
%qbshMatch: Match a given pitch vector for QBSH [對一個使用者輸入的音高向量進行 QBSH]
%	Usage: [distVec, scoreVec] = qbshMatch(songData, pitch, mrParam)
%		songData: song collection
%		pitch: input pitch vector
%		mrParam: MR parameters
%
%	The comparison results are written back into songData's field "score".

pitch=pitch(:)';				% Change to a row vector [改為列向量]
pitch=pvRestHandle(pitch, mrParam.useRest);		% Handle rests [處理休止符]
pitch=pitch-mean(pitch);			% Shift the pitch to have zero mean [平移到平均值為零]
pitchLen = length(pitch);			% Length of the pitch vector [音高向量長度（點數）]
maxSongLen = round(mrParam.lengthRatio*pitchLen);	% The song length is equal to mrParam.lengthRatio times the pitch length [標準歌曲長度最多只有輸入音高向量長度的 mrParam.lengthRatio 倍]
pitch=pitch(1:mrParam.pvrr:end);			% Down-sample to reduce computation [降點以減低運算量]
songNum=length(songData);
distVec=zeros(1, songNum);

switch mrParam.method
	case 'dtw1'
		for i=1:length(songData),
%			fprintf('%d/%d\n', i, length(songData));
			song=songData(i).pv;
			% ====== Take a reference song of an appropriate length [擷取歌曲適當長度]
			song=song(1:min(maxSongLen, length(song)));
			songLen = length(song);
			% ====== Shift song to have zero mean [平移歌曲，盡量使其和輸入音高向量有相同的基準（在此為平移至平均值為零）]
			songMean=mean(song(1:min(pitchLen, songLen)));	% Get mean of the song with the same length as the pitch vector [抓出歌曲的平均值（取其長度和輸入音高向量相同）]
			song=song-songMean;				% Shift to mean [平移到平均值]
			% ====== Down-sample to reduce computation [降點以減低運算量]
			song=song(1:mrParam.pvrr:end);
			% ====== Set DTW parameters [設定 DTW 參數]
			beginCorner=mrParam.beginCorner;	% 1: anchored beginning [頭固定]
			endCorner=mrParam.endCorner;		% 0: free end [尾浮動]
			% ====== Compute DTW [計算 DTW]
			distVec(i)=dtw4qbsh(pitch, song, mrParam);
			% ====== Compute DTW [計算 DTW]
		%	distVec(i) = dtw2mex(pitch, song, beginCorner, endCorner);
		end
	case 'dtw2'
		for i=1:length(songData),
%			fprintf('%d/%d\n', i, length(songData));
			song=songData(i).pv;
			% ====== Take a reference song of an appropriate length [擷取歌曲適當長度]
			song=song(1:min(maxSongLen, length(song)));
			songLen = length(song);
			% ====== Shift song to have zero mean [平移歌曲，盡量使其和輸入音高向量有相同的基準（在此為平移至平均值為零）]
			songMean=mean(song(1:min(pitchLen, songLen)));	% Get mean of the song with the same length as the pitch vector [抓出歌曲的平均值（取其長度和輸入音高向量相同）]
			song=song-songMean;				% Shift to mean [平移到平均值]
			% ====== Down-sample to reduce computation [降點以減低運算量]
			song=song(1:mrParam.pvrr:end);
			% ====== Set up DTW parameters [設定 DTW 參數]
			beginCorner=mrParam.beginCorner;	% 1: anchored beginning [頭固定]
			endCorner=mrParam.endCorner;		% 0: free end [尾浮動]
			% ====== Compute DTW (計算 DTW)
			distVec(i) = dtw2mex(pitch, song, beginCorner, endCorner);
		end
	case 'dtw1fixedPoint'
		for i=1:length(songData),
%			fprintf('%d/%d\n', i, length(songData));
			song=songData(i).pv;
			% ====== Take a reference song of an appropriate length [擷取歌曲適當長度]
			song=song(1:min(maxSongLen, length(song)));
			song=song(1:mrParam.pvrr:end);	% Down-sample to reduce computation [降點以減低運算量]
			% ====== Set up parameters for dtwFixedPoint
			param.dtwFunc='dtw1mex';
			param.anchorBeginning=mrParam.beginCorner;	% 1: anchored beginning [頭固定]
			param.anchorEnd=mrParam.endCorner;		% 0: free end [尾浮動]
			param.maxIterationNum=mrParam.dtwCount;
			[minDist, allDist, allPitchShift]=dtwFixedPoint(pitch, song, param);
			distVec(i)=minDist;
		end
	case 'dtw2fixedPoint'
		for i=1:length(songData),
%			fprintf('%d/%d\n', i, length(songData));
			song=songData(i).pv;
			% ====== Take a reference song of an appropriate length [擷取歌曲適當長度]
			song=song(1:min(maxSongLen, length(song)));
			song=song(1:mrParam.pvrr:end);	% Down-sample to reduce computation [降點以減低運算量]
			% ====== Set up parameters for dtwFixedPoint
			param.dtwFunc='dtw2mex';
			param.anchorBeginning=mrParam.beginCorner;	% 1: anchored beginning [頭固定]
			param.anchorEnd=mrParam.endCorner;		% 0: free end [尾浮動]
			param.maxIterationNum=mrParam.dtwCount;
			[minDist, allDist, allPitchShift]=dtwFixedPoint(pitch, song, param);
			distVec(i)=minDist;
		end	
	case 'ls'
		% ====== Set up LS parameters [設定 LS 參數]
		lowerRatio=mrParam.lowerRatio;	% 0.5;
		upperRatio=mrParam.upperRatio;	% 2.0;
		resolution=mrParam.resolution;	% 11;
		distanceType=1;
		[scaledVecSet, scaledVecLen]=scaledVecCreate(pitch, lowerRatio, upperRatio, resolution);
		for i=1:length(songData)
%			fprintf('%d/%d\n', i, length(songData));
			song=songData(i).pv;
			% ====== Take a reference song of an appropriate length [擷取歌曲適當長度]
			song=song(1:min(maxSongLen, length(song)));
			songLen = length(song);
			% ====== Down-sample to reduce computation (降點以減低運算量)
			song=song(1:mrParam.pvrr:end);
		%	% ====== Compute LS (計算 LS)
		%	distVec(i)=linScalingMex(pitch, song, lowerRatio, upperRatio, resolution, distanceType);
			% ====== Compute LS (計算 LS)
			distVec(i)=linScaling2Mex(scaledVecSet, scaledVecLen, song, distanceType);
		end
	case 'noteLevelDtw2'
		for i=1:length(songData),
%			fprintf('%d/%d\n', i, length(songData));
			% ====== Shift song to have zero mean [平移歌曲，盡量使其和輸入音高向量有相同的基準（在此為平移至平均值為零）]
			song=songData(i).pv;
			songLen = length(song);
			songMean=mean(song(1:min(pitchLen, songLen)));	% Get mean of the song with the same length as the pitch vector [抓出歌曲的平均值（取其長度和輸入音高向量相同）]
			songNote=songData(i).track(1:2:end)';		% Use pitch only [只取音高]
			songNote=songNote-songMean;		
			% ====== Note segmentation [進行音符切割]
			pitchTh=0.8;
			minNoteDuration=0.1;
			userNote = noteSegment(pitch, mrParam.timeStep, pitchTh, minNoteDuration);
			% ====== Use dtw2 to compute the distance between song and userNote. [使用 DTW2 來計算 song 和 userNote 之間的距離]
			userNote=userNote(1:2:end);		% Use pitch only [只取音高]
			beginCorner=mrParam.beginCorner;		% 1: anchored beginning [頭固定]
			endCorner=mrParam.endCorner;			% 0: free end [尾浮動]
			distVec(i) = dtw2mex(userNote, songNote, beginCorner, endCorner);
		end
	case 'dtw3'
		for i=1:length(songData),
%			fprintf('%d/%d\n', i, length(songData));
			% ====== Shift song to have zero mean [平移歌曲，盡量使其和輸入音高向量有相同的基準（在此為平移至平均值為零）]
			song=songData(i).pv;
			songLen = length(song);
			songMean=mean(song(1:min(pitchLen, songLen)));	% Get mean of the song with the same length as the pitch vector [抓出歌曲的平均值（取其長度和輸入音高向量相同）]
			songNote=songData(i).track(1:2:end)';		% Use pitch only [只取音高]
			songNote=songNote-songMean;
			% ====== Down-sample to reduce computation [降點以減低運算量]
			song=song(1:mrParam.pvrr:end);
			% ====== Set up DTW parametrs [設定 DTW 參數]
			beginCorner=mrParam.beginCorner;	% 1: anchored beginning [頭固定]
			endCorner=mrParam.endCorner;		% 0: free end [尾浮動]
			% ====== 計算 DTW
			distVec(i) = dtw3mex(pitch, songNote, beginCorner, endCorner);
		end
	otherwise
		error('Unknown method!');
end

scoreVec=100./(1+distVec);
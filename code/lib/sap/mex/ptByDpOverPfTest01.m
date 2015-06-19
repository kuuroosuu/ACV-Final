clear all; addMyPath;
fprintf('Compiling ptByDpOverPfMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility ptByDpOverPfMex.cpp /users/jang/c/lib/audio/audio.cpp /users/jang/c/lib/utility/utility.cpp
fprintf('Finish compiling ptByDpOverPfMex.cpp\n');
dos('copy /y ptByDpOverPfMex.mex* ..');
addpath ..
fprintf('Test ptByDpOverPfMex...\n');

waveFile='lately2.wav';		% 16 kHz, 16 bits
waveFile='../西路蟬聲唱.wav';
waveFile='但使龍城飛將在.wav';
%waveFile='yi_cuen_xiang_s_yi_cuen_huei.wav';
%waveFile='what_movies_have_you_seen_recently.wav';
%waveFile='what_would_you_like_to_know.wav';
%waveFile='/dataset/tangPoem/2007-音訊處理與辨識/92502067林子喬/二月黃鸝飛上林.wav';
waveFile='yankee_doodle.wav';		% 8K, 8bit
waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\霈宗_1\國家_張帝_1.wav';
waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\霈宗_1\不只是朋友_黃小琥_2.wav';
waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\霈宗_1\春夏秋冬_蔡小虎_1.wav';
waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\霈宗_1\堆積情感_黎明_1.wav';
%waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\Demo_0\好心情_李玟.wav';
%waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\Demo_1\兩隻老虎_不詳_0_s_2.wav';
%waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\renee_0\一見鍾情_藍心湄_1.wav';
%waveFile='yankee_doodle.wav';		% 8K, 8bit
waveFile='/dataSet/childSong/2007-音訊處理與辨識/19461108任佳王民/十個印第安人_不詳_0.wav';
waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\霈宗_1/冰雨_劉德華_1.wav';
waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave/吳明儒_1/月亮惹的禍_張宇_1_s.wav';
waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave/黃韋中_1/千里之外_周杰倫.費玉清_1_s.wav';
%waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\霈宗_1\一支小雨傘_洪榮宏_5.wav';
waveFile='D:\dataSet\dean\201005-mirRecording\Record_Wave\霈宗_1/冰雨_劉德華_1.wav';
wObj=waveFile2obj(waveFile);
if wObj.fs==16000
	wObj=wObj2wObj(wObj, 8000, 16, 1);
end

pfType=1;	% 0 for AMDF, 1 for ACF
ptOpt=ptOptSet(wObj.fs, wObj.nbits, pfType);
tic
%[pitch1, clarity, pfMat]=pitchTracking(wObj, ptOpt, 1);
[pitch1, clarity, pfMat]=pitchTracking(wObj, ptOpt, 1);
toc

frameNum=floor((length(wObj.signal)-ptOpt.overlap)/(ptOpt.frameSize-ptOpt.overlap));
frameTime=frame2sampleIndex(1:frameNum, ptOpt.frameSize, ptOpt.overlap)/wObj.fs;

ptOpt.useVolThreshold=0;
pitch2=pitchTracking(wObj, ptOpt);

ptOpt.useClarityThreshold=0;
pitch3=pitchTracking(wObj, ptOpt);

subplot(5,1,3);
h=get(gca, 'child'); delete(h(1));
pitch1(pitch1==0)=nan;
pitch2(pitch2==0)=nan;
pitch3(pitch3==0)=nan;
line(frameTime, pitch3-6, 'marker', '.', 'color', 'r');
line(frameTime, pitch2-3, 'marker', '.', 'color', 'g');
line(frameTime, pitch1,  'marker', '.', 'color', 'b');
title('Unbroken pitch (red), with clarity thresholding (green), with volume & clarity thresholding (blue)');
axis([-inf inf 30 84]);

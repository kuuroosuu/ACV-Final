clear all; addMyPath;
fprintf('Compiling ptByDpOverPfMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility ptByDpOverPfMex.cpp /users/jang/c/lib/audio/audio.cpp /users/jang/c/lib/utility/utility.cpp
fprintf('Finish compiling ptByDpOverPfMex.cpp\n');
%dos('copy /y ptByDpOverPfMex.mex* ..');
%addpath ..
fprintf('Test ptByDpOverPfMex...\n');

waveFile='yankee_doodle.wav';
waveFile='twinkle_twinkle_little_star.wav';
%waveFile='10LittleIndians.wav';	% Noisy!
wObj=waveFile2obj(waveFile);
pvFile=[wObj.file(1:end-3), 'pv']; if exist(pvFile), wObj.tPitch=asciiRead(pvFile); end
ptOpt=ptOptSet(wObj.fs, wObj.nbits);
pitch=pitchTracking(wObj, ptOpt, 1);
return

if wObj.amplitudeNormalized
	wObj.signal=wObj.signal*2^(wObj.nbits-1);
	wObj.amplitudeNormalized=0;
end
[pitch, clarity, pfMat, dpPath]=ptByDpOverPfMex(wObj.signal, wObj.fs, wObj.nbits, ptOpt);
clear all; addMyPath;
fprintf('Compiling ptByDpOverPfMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility -I/users/jang/c/lib/wave ptByDpOverPfMex.cpp /users/jang/c/lib/audio/audio.cpp /users/jang/c/lib/utility/utility.cpp /users/jang/c/lib/wave/waveRead4pda.cpp -output ptByDpOverPfMex
fprintf('Finish compiling ptByDpOverPfMex.cpp\n');
dos('copy /y ptByDpOverPfMex.mex* ..');
addpath ..
fprintf('Test ptByDpOverPfMex...\n');

load ttt.mat
[pitch, clarity, pfMat, dpPath]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
plot(pitch);

return

for i=1:1000
	fprintf('%d/%d\n', i, 1000);
	[pitch, clarity, pfMat, dpPath]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
	plot(pitch); drawnow
end

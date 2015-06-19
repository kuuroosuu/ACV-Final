% frameFlip.m 和 frameFlipMex.dll 的同步測試

fprintf('Compiling frameFlipMex.cpp...\n');
mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave frameFlipMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp -output frameFlipMex.dll

n=10000;
for i=1:n
	fprintf('%d/%d\n', i, n);
	frame=round(100*randn(256, 1));
	frame1=frameFlip(frame);
	frame2=frameFlipMex(frame);
	if sum(abs(frame1-frame2))~=0
		fprintf('Error!\n');
		plot([frame1, frame2]);
		return;
	end
end
addpath d:/users/jang/matlab/toolbox/sap

clear all
waveFiles={'yankeeDoodle_8k8b.wav', 'ivegotrhythm_16k16b.wav'};
waveFiles={'七月七日長生殿.wav'};	% This file causes frame2acfIntMex to overflow

% ====== Test AMDF
fprintf('AMDF test:\n');
for k=1:length(waveFiles)
	waveFile=waveFiles{k};
	[y, fs, nbits]=wavReadInt(waveFile);
	frameSize=256*fs/8000; overlap=0;
	frameMat=enframe(y, frameSize, overlap); 
	frameNum=size(frameMat, 2);
	fprintf('\twaveFile=%s, fs=%d, nbits=%d, frameSize=%d, overlap=%d, frameNum=%d\n', waveFile, fs, nbits, frameSize, overlap, frameNum);
	for method=1:3
		for i=1:frameNum
			frame=frameMat(:, i);
			out1=frame2amdf(frame, frameSize/2, method);
			out2=frame2amdfMex(frame, frameSize/2, method);
			difference(i)=max(abs(out1-out2));
		end
		fprintf('\t\tMethod=%d, max difference=%g\n', method, max(difference));		
	end
	fprintf('\t(The error is due the the fixed-point division of AMDF.)\n');
end

% ====== Test ACF
fprintf('ACF test:\n');
for k=1:length(waveFiles)
	waveFile=waveFiles{k};
	[y, fs, nbits]=wavReadInt(waveFile);
	frameSize=256*fs/8000; overlap=0;
	frameMat=enframe(y, frameSize, overlap); 
	frameNum=size(frameMat, 2);
	fprintf('\twaveFile=%s, fs=%d, nbits=%d, frameSize=%d, overlap=%d, frameNum=%d\n', waveFile, fs, nbits, frameSize, overlap, frameNum);
	for method=1:3
		for i=1:frameNum
			frame=frameMat(:, i);
			out1=frame2acf(frame, frameSize/2, method);
			out2=frame2acfMex(frame, frameSize/2, method);
			difference(i)=max(abs(out1-out2));
		end
		fprintf('\t\tMethod=%d, max difference=%g\n', method, max(difference));
	end
end

% ====== Test NSDF
fprintf('NSDF test:\n');
for k=1:length(waveFiles)
	waveFile=waveFiles{k};
	[y, fs, nbits]=wavReadInt(waveFile);
	frameSize=256*fs/8000; overlap=0;
	frameMat=enframe(y, frameSize, overlap); 
	frameNum=size(frameMat, 2);
	fprintf('\twaveFile=%s, fs=%d, nbits=%d, frameSize=%d, overlap=%d, frameNum=%d\n', waveFile, fs, nbits, frameSize, overlap, frameNum);
	for method=1:3
		for i=1:frameNum
			frame=frameMat(:, i);
			out1=frame2nsdf(frame, frameSize/2, method);
			out2=frame2nsdfMex(frame, frameSize/2, method);
			difference(i)=max(abs(out1-out2));
		end
		fprintf('\t\tMethod=%d, max difference=%g\n', method, max(difference));
	end
end
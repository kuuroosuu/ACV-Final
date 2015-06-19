wavfile = 'rain2.wav';
[y, fs] = wavread(wavfile);

frameSize = 512;
overlap = 340;

frameMat1 = buffer2(y, frameSize, overlap);	% Create overlapping buffer matrix
frameMat2 = buffer2mex(y, frameSize, overlap);
fprintf('Max absolute deviation is %g\n', max(max(abs(frameMat1-frameMat2))));

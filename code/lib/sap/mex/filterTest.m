close all; clear all;
fprintf('Compiling filterMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility filterMex.cpp /users/jang/c/lib/audio/audio.cpp /users/jang/c/lib/utility/utility.cpp

wavfile = 'yankeeDoodle_8k8b.wav';
[y, fs] = wavread(wavfile);

%order = 5;
%b = fir1(order, 900/(fs/2));    % Hamming is used by default

% ====== Test when the input is a vector 
y1 = filter([1 -1], [1 -0.98], y);
y2 = filterMex([1 -1], [1 -0.98], y);
fprintf('Length of y is %g\n', length(y));
maxDiff=max(max(abs(y1-y2)));
fprintf('Max absolute deviation is %g\n', maxDiff);
subplot(211); plot(1:length(y1), y1, 1:length(y2), y2);
subplot(212); plot(1:length(y1), abs(y1-y2));

% ====== Test when the input is a matrix
frameMat=buffer(y, 256, 0);
y1 = filter([1 -1], [1 -0.98], frameMat);
y2 = filterMex([1 -1], [1 -0.98], frameMat);
fprintf('Length of y1 is %g\n', length(y(:)));
maxDiff=max(max(abs(y1-y2)));
fprintf('Max absolute deviation is %g\n', maxDiff);
figure;
subplot(211); plot(1:length(y1(:)), y1(:), 1:length(y2(:)), y2(:));
subplot(212); plot(1:length(y1(:)), abs(y1(:)-y2(:)));
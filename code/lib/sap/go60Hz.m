waveFile='60Hz_noise.wav';
[y, fs, nbits]=wavReadInt(waveFile);
frameSize=3200;
overlap=0;
frameMat=buffer2(y, frameSize, overlap);
frameNum=size(frameMat, 2);
frame=frameMat(:, 1);
frame=y;

%subplot(2,1,1);
%time=(1:length(y))/fs;
%plot(time, y);
%set(gca, 'xlim', [min(time), max(time)]);
%xlabel('Time (sec)'); title(sprintf('Waveform of "%s"', strrep(waveFile, '_', '\_')));

% ====== By least-squares
f=60;
t=(1:length(frame))'/fs;
A=[cos(2*pi*f*t), sin(2*pi*f*t)];
x=A\frame;
p=x(1); q=x(2);
a=sqrt(p^2+q^2);
theta=atan(q/p);
estimated=-a*cos(2*pi*f*t-theta);	% Why do I need to put a minus sign here??? The model is a*cos(2*pi*f*t-theta).
error=norm(frame-estimated)/sqrt(frameSize);
fprintf('error=%f\n', error);
subplot(2,2,1);
plot(t, frame, t, estimated);
subplot(2,2,3);
plot(t, frame-estimated);

% ====== By FFT
[magSpec, phaseSpec, freq, powerSpecInDb]=fftOneSide(frame, fs);
[maxValue, maxIndex]=max(magSpec);
peakFreq=freq(maxIndex);
Y=fft(frame);
Y(maxIndex)=0; Y(length(frame)-maxIndex+2)=0;
estimated=ifft(Y);
error=mean(abs(frame-estimated));
error=norm(frame-estimated)/sqrt(frameSize);
fprintf('error=%f\n', error);
subplot(2,2,2);
plot(t, frame, t, estimated);
subplot(2,2,4);
plot(t, frame-estimated);

return

figure;
[mag, phase, freq, powerSpec]=fftOneSide(y, fs, 1);

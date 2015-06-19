function display(w)
% display: Display a wave object

fprintf('%s.fileName=%s\n', inputname(1), w.fileName);
fprintf('%s.fs=%g\n', inputname(1), w.fs);
fprintf('%s.nbits=%g\n', inputname(1), w.nbits);
fprintf('%s.frameSize=%g\n', inputname(1), w.frameSize);
fprintf('%s.overlap=%g\n', inputname(1), w.overlap);

subplot(2,1,1);
time=(1:length(w.signal))/w.fs;
plot(time, w.signal);
waveMax=2^w.nbits/2;
waveMin=-waveMax;
axis([-inf inf waveMin, waveMax]);
title(w.fileName);
xlabel('Time (seconds)');
grid on
subplot(2,1,2);
image(w.frameMatrix);
axis xy
title(sprintf('frameSize=%g, overlap=%d', w.frameSize, w.overlap));
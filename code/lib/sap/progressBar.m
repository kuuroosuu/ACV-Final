% This is a demo showing how to play a sound while displaying the progress bar
% Roger Jang, 20100604

waveFile='但使龍城飛將在.wav';
wObj=waveFile2obj(waveFile);
time=(1:length(wObj.signal))/wObj.fs;
h=plot(time, wObj.signal);
axis([-inf inf -1 1]);

sound(wObj.signal, wObj.fs);

h=line(nan*[1 1], [-1 1], 'color', 'r', 'erase', 'xor', 'linewidth', 5);
delta=20;
tic
for i=1:delta:length(wObj.signal);
	set(h, 'xdata', i/wObj.fs*[1 1]);
	drawnow;
end
toc

fprintf('You need to adjust delta such that the progress time should be as close as possible to wave duration.\n');
fprintf('Wave duration: %f\n', length(wObj.signal)/wObj.fs);
fprintf('Progress time: %f\n', toc);
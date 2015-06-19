if ~exist('pitch')
	fprintf('You need to run the system first to have the pitch to play!\n');
else
	fprintf('Playing the pitch...\n');
	pvPlay(freq2pitch(pitch.signals.values(:)), 11025/256);
end
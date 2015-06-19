function note = pv2noteStrict(pv, frameRate)
% pv2noteStrict: PV to note transformation
%	Usage: note = pv2noteStrict(pv, frameRate)
%		pv = [semitone1, semitone2, ...] at frameRate
%			(Default value of frameRate is 8000/256 second.)
%		note = [semitone1, time1, semitone2, time2, ...] with time
%			unit of 1 second

%	Roger Jang, 20001209, 20070531

if nargin<2, frameRate = 8000/256; end

prev = pv(1);
time = 1;
note = [];
for i=2:length(pv),
	if pv(i)==prev,
		time = time+1;
	else
		note = [note, prev, time];
		prev = pv(i);
		time = 1;
	end
end

note = [note, prev, time];
note(2:2:end) = note(2:2:end)/frameRate;
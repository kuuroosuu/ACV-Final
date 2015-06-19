function note = pv2note(pv, opt, showPlot)
% pv2note: Note segmentation on a given pitch vector
%	Usage:
%		note = pv2note(pv, opt, showPlot)
%
%	Description:
%		note = pv2note(pv, opt, showPlot) returns the note after performing note segmentation over a given pitch vector.
%			pv: The given pitch vector
%			opt: Options for the function
%				opt.method: 'simple' for a simple method, 'dp' for method of dynamic programming
%				opt.frameRate: frame rate (or pitch rate)
%				opt.minPitchError: Minimum pitch error (min. pitch deviation for each pitch point) within a note
%				opt.minNoteDuration: Minimum of note duration
%			showPlot: 1 for plotting
%			note: The output note sequence of [pitch, duration, pitch, duration, ...], where pitch is in semitone and duration is in sec.
%
%	Example:
%		pv=[49.21841 48.33041 48.33041 47.76274 48.04425 48.04425 47.76274 48.33041 48.04425 48.04425 48.04425 47.48574 47.48574 47.48574 47.48574 47.48574 47.76274 48.04425 48.04425 47.76274 47.76274 47.76274 47.21309 47.21309 47.21309 47.21309 47.21309 51.48682 51.48682 51.48682 51.48682 51.83658 51.83658 51.83658 51.83658 51.83658 54.92247 54.09792 54.09792 54.09792 54.09792 54.09792 54.09792 54.92247 54.92247 54.92247 54.92247 54.50529 54.50529 54.92247 54.50529 54.50529 54.09792 54.50529 55.34996 54.92247 54.92247 54.50529 54.50529 54.50529 54.50529 54.50529 54.50529 54.50529 54.50529 54.50529 54.50529 53.31086 53.31086 53.31086 53.31086 56.23796 55.78827 55.78827 56.23796 56.69965 56.69965 56.69965 57.17399 56.69965 56.69965 56.69965 56.69965 56.69965 56.69965 56.69965 57.17399 57.17399 57.17399 56.69965 56.69965 56.69965 56.69965 56.69965 56.69965 56.69965 59.76274 59.76274 59.76274 59.21309 59.21309 59.21309 58.68037 58.68037 58.68037 58.68037 54.09792 54.09792 54.09792 54.50529 54.50529 54.09792 54.09792 54.50529 54.92247 54.50529 54.50529 54.09792 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 53.69992 50.47805 51.48682 51.83658 52.19354 52.55803 52.19354 52.55803 51.83658 51.83658 52.55803 52.93035 52.93035 52.93035 52.93035 52.93035 51.14399 51.14399 54.50529 53.31086 52.55803 52.19354 52.19354 52.19354 52.55803 52.93035 54.09792 54.50529 54.92247 55.78827 56.23796 56.23796 55.78827 55.34996 54.09792 54.09792 54.09792 51.48682 50.15444 50.15444 50.80782 51.14399 51.14399 51.14399 51.14399 52.19354 52.19354 51.83658 51.83658 51.83658 51.48682 51.48682 51.48682 51.83658 51.83658 51.48682 51.48682 51.48682 51.48682 51.48682 50.80782 50.80782 52.55803 51.48682 51.14399 50.80782 51.14399 51.48682 51.48682 51.48682 50.80782 50.80782 50.80782 48.91732 48.62138 48.33041 48.33041 48.33041 48.04425 48.91732 48.91732 48.91732 49.21841 49.21841 48.91732 48.62138 48.33041 48.33041 48.33041 49.83678 48.62138 48.62138 48.62138 48.62138 48.62138 48.91732 49.52484 49.52484 48.91732 48.62138 48.33041];
%		waveFile='twinkle_twinkle_little_star.wav';
%		wObj=waveFile2obj(waveFile);
%	%	fprintf('Playback of the wave file...\n');
%	%	sound(wObj.signal, wObj.fs);
%		ptOpt=ptOptSet(wObj.fs, wObj.nbits);
%		pv=pitchTrack(wObj, ptOpt);
%	%	segment=segmentFind(pv); pv=pv(segment(2).begin:segment(2).end);
%		opt=pv2note('defaultOpt');
%		opt.frameRate=wObj.fs/(ptOpt.frameSize-ptOpt.overlap);
%		opt.method='simple';
%		subplot(211); note=pv2note(pv, opt, 1);
%		opt.method='dp';
%		subplot(212); note=pv2note(pv, opt, 1);
%	%	fprintf('Playback of the given pitch vector...\n');
%	%	pvPlay(pv, opt.frameRate);
%	%	fprintf('Playback of the identified notes...\n');
%	%	notePlay(note, 1, 1); 

%	Roger Jang, 20060527, 20070531, 20121029

if nargin<1, selfdemo; return; end
if ischar(pv) && strcmpi(pv, 'defaultOpt')    % Set the default options
	note.frameRate=8000/256;
	note.minPitchError=0.2;
	note.minNoteDuration=0.2;
	note.method='dp';
	return
end
if nargin<2||isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end

segment=segmentFind(pv);	% Segments are separated by silence

switch(lower(opt.method))
	case 'simple'
		for j=1:length(segment)
			thePv=pv(segment(j).begin:segment(j).end);
			note=[];
			pvInCurrentNote=thePv(1);		% PV in the current note
			k=1;
			for i=2:length(thePv)
				noteDuration=length(pvInCurrentNote)/opt.frameRate;
				nextPv=[pvInCurrentNote, thePv(i)];
				meanPitchError=mean(abs(nextPv-median(nextPv)));
				if meanPitchError>opt.minPitchError & noteDuration>opt.minNoteDuration		% Create a new note if the pitch diff. is high and note duration is long...
					note(k).pitch=median(pvInCurrentNote);
					note(k).duration=noteDuration;
					if k==1
						note(k).pvIndex=1:i-1;
					else
						note(k).pvIndex=note(k-1).pvIndex(end)+1:i-1;
					end
					pvInCurrentNote=thePv(i);
					k=k+1;
				else
					pvInCurrentNote=[pvInCurrentNote, thePv(i)];
				end
			end
			% === Process the last note
			note(k).pitch=median(pvInCurrentNote);
			note(k).duration=length(pvInCurrentNote)/opt.frameRate;
			if k==1
				note(k).pvIndex=1:length(thePv);
			else
				note(k).pvIndex=note(k-1).pvIndex(end)+1:length(thePv);
			end
			% === Merge the last two notes if the last note is too short
			if note(end).duration<opt.minNoteDuration
				if length(note)==1, error('A single note that is too short!'); end
				note(end-1).pvIndex=[note(end-1).pvIndex, note(end).pvIndex];
				note(end)=[];
				note(end).pitch=median(thePv(note(end).pvIndex));
				note(end).duration=length(note(end).pvIndex)/opt.frameRate;
			end
			% === Assign to segment
			segment(j).note=note;
		end
	case 'dp'
		for j=1:length(segment)
			thePv=pv(segment(j).begin:segment(j).end);
			opt2.segmentNum=6;
			opt2.minNoteDuration=0.1;
			opt2.frameRate=opt.frameRate;
			for i=1:length(thePv)
				opt2.segmentNum=i;
				[normalizedMinDist, dtwPath, dtwTable, note]=sequenceSegmentViaDp(thePv, opt2);
				if normalizedMinDist<=opt.minPitchError, break; end
			end
			segment(j).note=note;
		end
	otherwise
		error('Unknown method %s in %s!', opt.method, mfilename);
end

if showPlot
	pv2=pv; pv2(pv2==0)=nan;
	plot((1:length(pv2))/opt.frameRate, pv2, '.-g');
	hold on;
	for i=1:length(segment)
		notePlot(segment(i).note, 1, 'b', segment(i).begin/opt.frameRate);
	end
	hold off
	title(sprintf('Method=%s, no. of notes=%d', opt.method, length([segment.note])));
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);

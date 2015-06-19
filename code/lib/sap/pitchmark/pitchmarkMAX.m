function pitchmarks = pitchmarkMAX(y, fs, maxPeriod, minPeriod, SFS_pitchcurve, plotopt)

%PITCHMARKMAX Pitch Mark
%   Usage : pitchmarks = pitchmarkMAX(y, fs, maxPeriod, minPeriod, SFS_pitchcurve, plotopt)
%   Selfdemo : pitchmarkMAX;

if nargin <6, plotopt = 0; end;
if nargin <3, maxPeriod = 0.007; minPeriod = 0.0065; SFS_pitchcurve = 200; plotopt = 1; end;
if nargin==0, selfdemo; return; end;

pitch = 1/((maxPeriod+minPeriod)/2); %︳衡ㄓpitch

% ===== т程y
[maximum, pitchmarks] = max(y);
% ┕тN pitch marks
while (1),
   currentIndex = pitchmarks(end);
   minIndex = currentIndex+floor(minPeriod*fs);
   maxIndex = currentIndex+ceil(maxPeriod*fs);
   if maxIndex >= length(y) ,
      if (length(y)- pitchmarks(end)) > (fs/pitch), %安程pitch markび
         begin = pitchmarks(end)+ round((fs/pitch)/2);
         candidateY = y(begin : end);
         [junk, newIndex] = max(candidateY);
         pitchmarks = [pitchmarks, begin+newIndex-1];
      end;
      break; 
   end;
   candidateY = y(minIndex:maxIndex);
   [junk, newIndex] = max(candidateY);
   pitchmarks = [pitchmarks, minIndex+newIndex-1];
end
%┕玡тM pitch marks
while (1)
   currentIndex = pitchmarks(1);
   maxIndex = currentIndex-floor(minPeriod*fs);
   minIndex = currentIndex-ceil(maxPeriod*fs);
   if minIndex <= 0 , %临Τ程pitch mark
      if pitchmarks(1) > (fs/pitch), %安ヘ玡材pitch mark
         tail = pitchmarks(1) - round((fs/pitch)/2);
         candidateY = y(1:tail);
         [junk, newIndex] = max(candidateY);
         pitchmarks = [newIndex, pitchmarks];
      end;
      break;
   end;
   candidateY = y(minIndex:maxIndex);
   [junk, newIndex] = max(candidateY);
   pitchmarks = [minIndex+newIndex-1, pitchmarks];
end;

%タpitch mark岿粇
pitchmarks((diff(pitchmarks)==0)) = [];

if plotopt,
   plotpitchmarks(y, fs, pitchmarks, SFS_pitchcurve);
end;


% ============================================================== %
% selfdemo
% ============================================================== %
function selfdemo
%wavefile = 'test.wav';
wavefile = 'C:\Users\tiny\Desktop\tiny\dataset\target\cut\バ_4.wav';
plotopt = 1;
[y, fs] = wavread(wavefile);

% 璸衡 pitch
SFS_pitchcurve = wave2pitch_SFS(wavefile);
SFS_pitchcurve(find(SFS_pitchcurve<50 | SFS_pitchcurve>600)) = [];
if isempty(SFS_pitchcurve), error('Pitch is not available!'); end;

% 倒﹚ pitchmarkMIN 把计
PitchMean = mean(SFS_pitchcurve);
a = (1/PitchMean);
b = a*0.5;
maxPeriod = a + b;
minPeriod = a - b;

pitchmarks = feval(mfilename, y, fs, maxPeriod, minPeriod, SFS_pitchcurve, plotopt);

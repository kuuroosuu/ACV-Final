function  pmarkPrm = getpmarkPrm(y, fs, plotopt, SFS_pitchcurve)

% GETPMARKPRM Get parameters for pitch marking function.
%   Usage : pmarkPrm = getpmarkPrm(y, fs, plotopt, SFS_pitchcurve)

if nargin == 0, 
   [y, fs] = wavread('la.wav');
   plotopt = 1; 
   SFS_pitchcurve = [176   163   150   147   146   145   146   150   144   145   144   144   143   142   142   142   142   142,...
         142   141   142   141   140   140   140   140   141   140   141   142   141   141   141   141   143   144,...
         146   146   146   149   149   157   155   155   147   167]; 
end;
if nargin < 4, 
   wavwrite(y, fs, 16, 'temp.wav');
   SFS_pitchcurve = wave2pitch_SFS([pwd '\temp.wav']);
   SFS_pitchcurve((SFS_pitchcurve<50 | SFS_pitchcurve>600)) = [];
   if isempty(SFS_pitchcurve), error('Pitch is not available!'); end;   
end;
if nargin < 3, plotopt = 0; end;

PitchMean = mean(SFS_pitchcurve);
a = (1/PitchMean);
b = a*0.7;
maxPeriod = a + b;
minPeriod = a - b;
PeakMarks = pitchmarkMAX(y, fs, maxPeriod, minPeriod, SFS_pitchcurve, plotopt); %MAX
ValleyMarks = pitchmarkMIN(y, fs, maxPeriod, minPeriod, SFS_pitchcurve, plotopt); %MIN
PeakMarks_pitchcurve = fs./diff(PeakMarks)'; %row-based vector.
ValleyMarks_pitchcurve = fs./diff(ValleyMarks)'; %row-based vector.
SFS_pitchcurve1 = scalevector(SFS_pitchcurve, length(PeakMarks_pitchcurve))'; %row-based vector
SFS_pitchcurve2 = scalevector(SFS_pitchcurve, length(ValleyMarks_pitchcurve))'; %row-based vector
dpA = abs(PeakMarks_pitchcurve - SFS_pitchcurve1);
dpB = abs(ValleyMarks_pitchcurve - SFS_pitchcurve2);
A = sum(dpA)/length(PeakMarks_pitchcurve);
B = sum(dpB)/length(ValleyMarks_pitchcurve);
pmarkPrm.meanPitch = mean(SFS_pitchcurve);
pmarkPrm.maxPeriod = maxPeriod;
pmarkPrm.minPeriod = minPeriod;
pmarkPrm.candNum   = 3; %Candidates Number

if A < B, %Select MAX is better.
   pmarkPrm.PorV = 0;
   pmarkPrm.pitchcurve = [SFS_pitchcurve1; SFS_pitchcurve1(end)]; %Add additional element for the using of DP.
else      %Select MIN is better.
   pmarkPrm.PorV = 1;
   pmarkPrm.pitchcurve = [SFS_pitchcurve2; SFS_pitchcurve2(end)]; %Add additional element for the using of DP.
end;
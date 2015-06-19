function pitchmarks = pitchmarkdp(y, fs, plotopt, pmarkPrm)
% PITCHMARKDP Pitch marking based on Dynamic Programming.
%   Usage : pitchmarks = pitchmarkdp(y, fs, plotopt, pmarkPrm)
%
%   pmarkPrm is a struct variable.
%   You can type 'pmarkPrm = getpmarkPrm(y, fs);' to get this variable.
%   or use 'pmarkPrm = getpmarkPrm(y, fs, plotopt, SFS_pitchcurve);' to get more precise parameters.
%
%   Selfedemo : pitchmarkdp
%
% Cheng-Yuan Lin 2004, April, 5

if nargin == 0, selfdemo; return; end;
if nargin <4, pmarkPrm = getpmarkPrm(y, fs); end;

%Setting Parameters.
meanPitch = pmarkPrm.meanPitch;
maxPeriod = pmarkPrm.maxPeriod;
minPeriod = pmarkPrm.minPeriod;
candNum = pmarkPrm.candNum;
PorV = pmarkPrm.PorV;
SFS_pitchcurve = pmarkPrm.pitchcurve;

%Main Program.
sampleIndex = getlocation(y, fs, meanPitch, candNum, PorV, maxPeriod, minPeriod, 0);
probmatrix  = getstateprob(y, sampleIndex);
pitchmarks  = rundp(probmatrix, sampleIndex, SFS_pitchcurve, fs);

if plotopt,   
   plotpitchmarks(y, fs, pitchmarks, SFS_pitchcurve);
end;

% ============================================================ %
%                Sub function: getstateprob                    %
% ============================================================ %
function probmatrix = getstateprob(y, sampleIndex)
row = size(sampleIndex, 1);
col = size(sampleIndex, 2);
probmatrix = ones(row, col);
dist = zeros(1, row);  %dist is from 0 ~ 2, and prob is from 1 ~ 0.
scale = 1;
for c = 1 : col,
   dist(2:3) = y(sampleIndex(2:end, c)) - y(sampleIndex(1, c));
   dist(2:3) = dist(2:3)*scale;
   prob = 1-abs(dist/2);
   probmatrix(:, c) = log10(prob)';
   %probmatrix(:, c) = prob';
end;

% ============================================================ %
%                   Sub function: rundp                        %
% ============================================================ %
function pmarkdp = rundp(probmatrix, sampleIndex, SFS_pitchcurve, fs)
row = size(sampleIndex, 1);
col = size(sampleIndex, 2);
Indexmatrix = ones(row, col);
maxDist = mean(SFS_pitchcurve);

%差一個pitch points
if (length(SFS_pitchcurve) == col-1 ), SFS_pitchcurve = [SFS_pitchcurve SFS_pitchcurve(end)]; end;
if col > length(SFS_pitchcurve), col = length(SFS_pitchcurve); end;

for c = 1 : col-1,
   for r = 1 : row,
      temp = sampleIndex(r, c+1) - sampleIndex(:, c);
      pitch = fs./temp;
      %根據pitch curve算出與目前pitch的機率值
      %logprob = -(abs(temp-SFS_pitchcurve(c+1))./temp);

      dist = abs(pitch - SFS_pitchcurve(c+1));
      dist(find(dist>maxDist)) = maxDist;
      prob = -(dist/maxDist)+1;
      logprob = log10(prob)';
      %logprob = prob';

      logprob = logprob + probmatrix(:, c)';
      [cumvalue, index] = max(logprob);

      probmatrix(r, c+1) = probmatrix(r, c+1) + cumvalue;
      Indexmatrix(r, c+1) = index;
   end;
end;

%Back tracking.
backIndex = ones(1, col); %backIndex紀錄的是前一個index的值, 例如Indexmatrix 第17行紀錄的是 第16個pamrk的index
pmarkdp   = ones(1, col);
[junk, backIndex(end)] = max(probmatrix(:,end));
pmarkdp(end) = sampleIndex(backIndex(end), end);
for k = col-1 : -1 : 1,
   pmarkdp(k) = sampleIndex(backIndex(k+1), k);
   backIndex(k) = Indexmatrix(backIndex(k+1), k);
end;

% ============================================================== %
% selfdemo
% ============================================================== %
function selfdemo
[y, fs] = wavread('la.wav'); y = y(4000:7400);
%[y, fs] = wavread('E:\ProgramCode\DataSet\TTS_wavedata\CSWaveFiles\11021.wav');
plotopt = 1;
pmarkPrm = getpmarkPrm(y, fs, plotopt);
pitchmarks = feval(mfilename, y, fs, plotopt, pmarkPrm);
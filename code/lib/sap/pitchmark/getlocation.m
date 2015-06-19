function sampleIndex = getlocation(y, fs, pitch, candNum, PorV, maxPeriod, minPeriod, plotopt)

%GETLOCATION Find peaks or valleys location to search best pitch marks.
%   Usage : sampleIndex = getlocation(y, fs, pitch, candNum, PorV, maxPeriod, minPeriod, plotopt)
%   Selfdemo : getlocation;
%   
% Cheng-Yuan Lin

if nargin <5, plotopt = 0; end;
if nargin <3, maxPeriod = 0.007; minPeriod = 0.0065; plotopt = 1; end;
if nargin==0, selfdemo; return; end;

check = 0;
Xdistscale    = 20;  %籔程翴x禯瞒 +- 20翴钡
Ydistscale    = 0.2; %籔程翴y禯瞒ぶ0.2钡, normalizedぇ把计

if PorV == 0, %0 for peaks
   [junk, pmark] = max(y);
else %1 for valleys
   [junk, pmark] = min(y);
end;

minYhead = pmark-50; if minYhead <= 0, minYhead = 1; end;
minYtail = pmark+50; if minYtail > length(y), minYtail = length(y); end;
sampleIndex   = getIndex(y(minYhead:minYtail), minYhead, minYtail, candNum, Xdistscale, Ydistscale, PorV);
singlePeriod  = (fs/pitch);

%┕т
while (1),
    %check
   currentIndex = sampleIndex(end, 1);
   minIndex = currentIndex+floor(minPeriod*fs);
   maxIndex = currentIndex+ceil(maxPeriod*fs);
   if maxIndex >= length(y) ,
      if (length(y)- sampleIndex(end, 1)) > singlePeriod, %安程pitch markび
         begin = sampleIndex(end, 1) + round(singlePeriod/2);
         candidateY = y(begin : end);
         [tempIndex,check] = getIndex(candidateY, begin, length(y), candNum, Xdistscale, Ydistscale, PorV);
         sampleIndex = [sampleIndex; tempIndex];
      end;
      break; 
   end;
   candidateY = y(minIndex:maxIndex);
   [tempIndex,check] = getIndex(candidateY, minIndex, maxIndex, candNum, Xdistscale, Ydistscale, PorV);
   sampleIndex = [sampleIndex; tempIndex];
   if (check == 1)
      break;
   end
end;

%┕玡т
while (1)
   currentIndex = sampleIndex(1, 1);
   maxIndex = currentIndex-floor(minPeriod*fs);
   minIndex = currentIndex-ceil(maxPeriod*fs);
   if minIndex <= 0 , %临Τ程pitch mark
      if sampleIndex(1, 1) > singlePeriod, %安ヘ玡材pitch mark
         tail = sampleIndex(1, 1) - round(singlePeriod/2);
         candidateY = y(1:tail);
         [tempIndex,check] = getIndex(candidateY, 1, tail, candNum, Xdistscale, Ydistscale, PorV);
         sampleIndex = [tempIndex; sampleIndex];
      end;
      break;
   end;
   candidateY = y(minIndex:maxIndex);
   [tempIndex,check] = getIndex(candidateY, minIndex, maxIndex, candNum, Xdistscale, Ydistscale, PorV);
   sampleIndex = [tempIndex; sampleIndex];
   if (check == 1)
      break;
   end
end;

sampleIndex = sampleIndex';

% ===========Plot figure
if plotopt,
   figure;
   pmark = sampleIndex(1, 1:end-1);
   subplot(2,1,1);plot((1:length(y)), y, '.-');
   %礶 pitch mark
   for i=1:length(pmark),
      line(pmark(i), y(pmark(i)), 'marker', 'o', 'color', 'r');
   end
   axis([-inf inf -1 1]);
   axis([0.1 0.8 -1 1]);
   ylabel('Amplitude');
   title('Wave form');
   axis tight;
   legend('wave','pitch marks');
   
   subplot(2,1,2);
   dpmark = diff(pmark);
   pitch = fs./dpmark;
   plot(pitch, 'ro-');
end;


function [tempIndex,check] = getIndex(candidateY, minIndex, maxIndex, candNum, Xdistscale, Ydistscale, PorV)
check = 0;
if PorV == 0,
   [value, idx] = localmax(candidateY);
else % PorV == 1,
   [value, idx] = localmin(candidateY);
end;

%璝栋, 玥тいみ翴蠢
if isempty(idx), 
   idx = round((minIndex + maxIndex)/2 - minIndex);
   value = candidateY(idx);
end;

%т程翴(pitch mark程Τぇ竚), 荷秖тcandNum才兵ン翴
dist = abs(value - value(1));
dist = dist/( max(candidateY) - min(candidateY)); %Normalized.
idx   = idx(find(dist < Ydistscale));
if (numel(idx)~=0)
    Index = idx(1);
else
    Index = idx;
end
range = [Index - Xdistscale : Index + Xdistscale];
idx = idx(find(ismember(idx, range)));

%安idx计筁(禬筁candNum┮﹚), 玥candNum
plen = length(idx);
if plen>candNum, plen = candNum; end;
tempIndex = minIndex+idx(1:plen)-1;
if length(tempIndex)< candNum,
   cnts = candNum - length(tempIndex);
   if (length(tempIndex)~=0)
   tempIndex = [tempIndex tempIndex(end)*ones(1, cnts)];
   %check = 0;
   else
       check = 1;
   end
end;
%end

function selfdemo
wavefile = 'la.wav';
plotopt = 0;
% ====== Read wavefile
[y,fs] = wavread(wavefile);
y = y-mean(y);
% ====== Endpoint detection
y = y(3900:7300);
%======== Pitch tracking
[p0, p1, p2, p3, p4] = wave2pitch(y,fs,0);
p2 = p2(find(p2~=0));
pitch = nanmedian(p2);
if isempty(p2) | isempty(pitch), 
   fprintf('The wavefile has some problems!\n');
   return;
end;
pitch = smtn2frq(pitch);
a = (1/pitch);
b = a*0.5;
maxPeriod = a + b;
minPeriod = a - b;
pitch = mean(pitch);
candNum = 3;
PorV    = 1;
sampleIndex = feval(mfilename, y, fs, pitch, candNum, PorV, maxPeriod, minPeriod, 1);

function [rawpitch, pitch, framedY, aveEng] = ptrackfcn(y, fs, plotopt, pitchparam)

% PTRACKFCN : Pitch Tracking Function
%    Usage [rawpitch, pitch, framedY, aveEng] = ptrackfcn(y, fs, plotopt, pitchparam).
% 
% Gavins 2001/10/15
% Gavins 2002/9/17 refined.

if nargin==0, selfdemo; return; end;
if nargin<3, plotopt =1 ; end;
if nargin<4, 
   pitchparam.frameSize = 320;
   pitchparam.overlap = 0;
   pitchparam.maxFreq = 600;
   pitchparam.lowpassOrder = 5;
   pitchparam.clipThred = 0.25;
   pitchparam.maxThred = 0.1;
   pitchparam.engthred = 0.2;
end;

% Shift mean to zero.
y = y-mean(y);
% Send to low-pass filter (Hamming window used by default)
maxFreq = pitchparam.maxFreq; %Assumed that frequency of human's speech is most 1K Hz.
b = fir1(pitchparam.lowpassOrder, maxFreq/(fs/2));
y = filter(b, 1, y);

% Frame blocking then computing average energy and zero crossing rate.
% Pitch tracking : center clipping, localmaxima, period extraction.
frameSize = pitchparam.frameSize;
overLap   = pitchparam.overlap;
framedY   = buffer2mex(y, frameSize, overLap); %function buffer in Signal toolbox
frameNum  = size(framedY,2);
seqCRLen  = frameSize*2-1;
sequenceCR= zeros(frameNum,seqCRLen);
clipThred = pitchparam.clipThred;
maxThred  = pitchparam.maxThred;
engThred  = pitchparam.engthred;

% Computing energy and using "min+(max-min)*0.2" rule to define threshold.
aveEng = getaveEnergy(framedY,frameSize);
%[a,b]  = sort(aveEng);
%c = fliplr(b);
%engThred = min(aveEng)+(mean(aveEng(c(1:round(length(c)*0.2)))) - min(aveEng)) *0.2;
engThred = engThred*mean(aveEng)*0.1;

rawpitch = zeros(frameNum,1);
for i = 1 : frameNum,
   sequenceCR(i,:) = atcorr(framedY(:,i));
   sequenceCR(i,:) = centclip(sequenceCR(i,:), clipThred);
   sequenceCR(i,:) = localmax(sequenceCR(i,:), maxThred);
   %Energy太小的音框不考慮求取pitch
   if aveEng(i) > engThred , 
      rawpitch(i) = freqExtract(sequenceCR(i,:), fs);
   else
      rawpitch(i) = 0;
   end;
end;
pitch    = pitchTrim(rawpitch);
semitone = frq2smtn(pitch);

% ====== Plot pitch figure.
if plotopt,
   time = (1:frameNum)*(frameSize-overLap)/fs;
   
   figure;
   subplot(5,1,1); 
   plot((1:length(y))/fs, y);axis tight;
   xlabel('time');ylabel('Amplitude');
   title(['(PCM Signal)']);
  
   subplot(5,1,2);
   plot(time, aveEng);axis tight;
   xlabel('time');ylabel('Average Energy');
   
   subplot(5,1,3);
   plot(time, rawpitch,'b-',time,rawpitch,'r*');axis tight;
   xlabel('time');ylabel('Raw Pitch');
   
   subplot(5,1,4);
   plot(time,pitch(1:length(time)),'b-',time,pitch(1:length(time)),'r*');
   axis([-inf inf 0 max(pitch)]);
   xlabel('time');ylabel('Pitch refined');
   
   subplot(5,1,5);
   plot(time,semitone(1:length(time)),'b-',time,semitone(1:length(time)),'r*');
   axis([-inf inf 0 max(semitone)]);
   xlabel('time');ylabel('Pitch to Semitone');
end;

function selfdemo
filename = 'test.wav';
plotopt = 1 ;
[y, fs] = wavread(filename);
[pitch, semitone] = ptrackfcn(y, fs, plotopt);

% ====================================================%
%                   Sub function                      %
% ====================================================%
% Short term average energy
function aveEng = getaveEnergy(matrix,vectorLen)
%aveEng = log(sum(power(matrix,2))/vectorLen);
aveEng = sum(power(matrix,2));

% Short term autocorrelation
function sequenceCR = atcorr(vector)
vectorLen  = length(vector);
vectorCOPY = vector; %vector duplicate
sequenceCR = zeros(1,vectorLen);
for i = 1 : vectorLen,
   sequenceCR(i) = sum(vector(i:vectorLen).*vectorCOPY(1:vectorLen-i+1));
end;
sequenceCR = [fliplr(sequenceCR(2:end)) sequenceCR]; %Sequence must be symmetric

% Short term center clipping
function sequenceCR = centclip(vector,clipThred)
vectorLen = length(vector);
sequenceCR = zeros(1,vectorLen);
clipvalue = max(abs(vector))*clipThred;
temp = abs(vector)-clipvalue;
index = find(temp>0);
sequenceCR(index) = temp(index).*sign(vector(index));

% Short term localmaximum
function sequenceCR = localmax(vector,maxThred)
vectorLen = length(vector);
sequenceCR = zeros(1,vectorLen);
maxvalue = max(abs(vector))*maxThred;
b = [vector 0];
c = [0 vector];
index1 = find([b - c]>0);
index2 = find([b - c]<=0);
b(index1) = 1;
b(index2) = 0;
maxIndex  = find(diff(b)<0);
Index  = maxIndex(find(vector(maxIndex)>maxvalue));
sequenceCR(1,Index) = 1;

% Short term frequency extraction
function rawpitch = freqExtract(vector,fs)
period = max(diff(find(vector~=0)));
if isempty(period), rawpitch = 0; return; end;
rawpitch = 1/(period/fs);

% Pitch fine tuning(trimming)
function trimedPitch = pitchTrim(rawpitch)
temp = buffer(rawpitch, 5, 2);
new  = median(temp,1);
trimedPitch = reshape([new'*ones(1,3)]',length(new)*3,1);
whichmin = min(length(trimedPitch), length(rawpitch));
trimedPitch = trimedPitch(1:whichmin);
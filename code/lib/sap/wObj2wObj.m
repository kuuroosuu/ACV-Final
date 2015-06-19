function wObj2=wObj2wObj(wObj, targetFs, targetNbits, targetChannelNum)
% wObj2wObj: Wave obj to wave obj conversion
%
%	Example:
%		waveFile='what_movies_have_you_seen_recently.wav';
%		wObj1=waveFile2obj(waveFile);
%		targetFs=11000;
%		targetNbits=8;
%		targetChannelNum=1;
%		wObj2=wObj2wObj(wObj1, targetFs, targetNbits, targetChannelNum);
%		timeSpan=[0.25, 0.28];
%		index1=round(timeSpan(1)*wObj1.fs):round(timeSpan(2)*wObj1.fs);
%		index2=round(timeSpan(1)*wObj2.fs):round(timeSpan(2)*wObj2.fs);
%		plot(index1/wObj1.fs, wObj1.signal(index1), '.-', index2/wObj2.fs, wObj2.signal(index2), '.-'); 

%	Roger Jang, 20100918

if nargin<1, selfdemo; return; end
if nargin<3, targetFs=8000; end
if nargin<4, targetNbits=8; end
if nargin<5, targetChannelNum=1; end

% ====== If specs are the same, copy input to output directly.
wObj2=wObj;
if wObj2.fs==targetFs & wObj2.nbits==targetNbits & size(wObj2.signal,2)==targetChannelNum
	return;
end

% ====== Resample if fs is different.
if wObj2.fs~=targetFs
	if size(wObj2.signal, 2)==1
		column1=resample(wObj2.signal(:,1), targetFs, wObj2.fs);
		wObj2.signal=column1;
	else
		column1=resample(wObj2.signal(:,1), targetFs, wObj2.fs);
		column2=resample(wObj2.signal(:,2), targetFs, wObj2.fs);
		wObj2.signal=[column1, column2];
	end
	wObj2.fs=targetFs;
end

% ====== Nbits conversion if nbits is different and useIntSignal
if wObj2.nbits~=targetNbits
	if ~wObj2.amplitudeNormalized
		wObj2.signal=wObj2.signal*2^(targetNbits-wObj2.Nbits-1);
	end
	wObj2.nbits=targetNbits;
end

% ====== Channel conversion if necessary
if size(wObj2.signal, 2)==2 & targetChannelNum==1
	wObj2.signal=wObj2.signal(:,1);
end
if size(wObj2.signal, 2)==1 & targetChannelNum==2
	wObj2.signal=[wObj2.signal, wObj2.signal];
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);

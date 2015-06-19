function audio2=audioConvert(audio, fs2, nbits2, channel2)
% audioConvert: Convert the format of a given audio object
%
%	Usage:
%		audio2=audioConvert(audio, fs2, nbits2, channel)
%
%	Description:
%		Conversion of a given audio to the specified format
%
%	Example:
%		audioFile='youAtLeast.mp3';
%		audio=myAudioRead(audioFile);
%		audio2=audioConvert(audio, 8000, 16, 1);

%	Roger Jang, 2004xxxx, 20141218

if nargin<1, selfdemo; return; end
if nargin<3, fs2=8000; end
if nargin<4, nbits2=8; end
if nargin<5, channel=1; end

if isstr(audio), audio=myAudioRead(audio); end

fs=audio.fs;
nbits=audio.nbits;
channel=size(audio.signal, 2);
audio2=audio; audio2.fs=fs2; audio2.nbits=nbits2;

% 把计妓钡肚
if (fs==fs2) & (nbits==nbits2) & (channel2==channel)
	return
end

% 把计ぃ秈︽锣传
if channel==2 & channel2==1		% Conversion from stereo to mono
	audio2.signal=mean(audio2.signal, 2);
end
if fs~=fs2
	temp=[];
	for i=1:size(audio2.signal, 2)
		temp(:,i)=resample(audio2.signal(:,i), fs2, fs);		% Resample to fs2
	end
	audio2.signal=temp;
end
if channel==1 & channel2==2		% Conversion from mono to stereo
	audio2.signal=[audio2.signal; audio2.signal];
end

% We did not consider nbits???

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
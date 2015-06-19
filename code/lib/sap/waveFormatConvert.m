function waveObj=waveFormatConvert(waveObj, targetFs, targetNbits, targetChannelNum)
% waveFormatConvert: Wave format conversion between wave objects
%	Usage: waveObj=waveFormatConvert(waveObj, targetFs, targetNbits, targetChannelNum)

%	Roger Jang, 20100410

if nargin<2, targetFs=8000; end
if nargin<3, targetNbits=8; end
if nargin<4, targetChannelNum=1; end

sourceChannelNum=size(waveObj.signal, 2);
if sourceChannelNum>targetChannelNum
	waveObj.signal=waveObj.signal(:,1);			% Take a single channel only
end
if sourceChannelNum<targetChannelNum
	waveObj.signal=[waveObj.signal, waveObj.signal];	% Duplicate to have 2 channels
end

if waveObj.fs~=targetFs
	waveObj.signal=resample(waveObj.signal, targetFs, waveObj.fs);
	waveObj.fs=targetFs;
end

if waveObj.nbits~=targetNbits
	waveObj.signal=waveObj.signal*2^(targetNbits-waveObj.nbits);		% bit-resolution conversion
	waveObj.nbits=targetNbits;
end

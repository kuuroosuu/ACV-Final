function ad2wavFile(adFile, wavFile)
%ad2wavFile Read an AD file and save it as a WAV file
%	Usage: ad2wavFile(adFile, wavFile)
%	
%	For example:
%		ad2wavFile('A01.AD', 'test.wav')

%	Roger Jang, 20090603

nbits=16;
useByteSwap=1;
fs=16000;

y=rawRead('A01.AD', nbits, useByteSwap);
wavwrite(y/32768, 16000, wavFile);
fprintf('Saved %s\n', wavFile);

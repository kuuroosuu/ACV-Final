function w=wave(fileName)
%WAVE Wave class constructor.
%	Usage: w = wave(fileName)
%		creates a polynomial object from the wav file fileName.

%	Roger Jang, 20060815

p.fileName=fileName;
[p.signal, p.fs, p.nbits]=wavread(fileName);
p.signal=p.signal*2^p.nbits/2;		% Convert to integer
p.frameSize=0;
p.overlap=0;
p.frameMatrix=[];
w=class(p, 'wave');
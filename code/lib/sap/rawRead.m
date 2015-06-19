function y=rawRead(rawFile, nbits, useByteSwap)
% rawRead: Read a binary sound file dumped from AP170, or from some other sources (such as the SAM format)
%	Usage: y=rawRead(rawFile, nbits, useByteSwap)
%		If nbits=8, the range of y is [0, 255];
%		If nbits=16, the range of y is [-32768, 32767]
%		If useByteSwap=1, do byte swapping for 16-bit sound data
%	
%	For example:
%	
%	y=rawRead('D:\dataSet\KeelePitch\Speech\f1nw0000.pes', 16);
%	sound(y/(2^16/2), 20000);
%
%	y=rawRead('C:\dataSet\fda_database\sb\sb001.sig', 16, 1);
%	sound(y/(2^16/2), 20000);

%	Roger Jang, 20051202, 20070212

if nargin<2, nbits=8; end
if nargin<3, useByteSwap=0; end

fid=fopen(rawFile, 'r');
switch nbits
	case 8
		y=fread(fid, 'uchar');
	case 16
		y=fread(fid, 'short');
		if useByteSwap						% 處理 byte swapping
			for i=1:length(y)
				bits=dec2bin(y(i), 16);
				bits=strrep(bits, '/', '1');		% 遇到負數，會出現'/'，將'/'改成'1'，就是2's complement
				bits=[bits(9:16), bits(1:8)];		% byte swapping
				signBit=bits(1);
				if signBit=='0';
					y(i)=bin2dec(bits);		% 正數
				else
					y(i)=-(65535-bin2dec(bits)+1);	% 負數
				end
			end
		end
	otherwise
		error('Unsupported number of bits!');
end

fclose(fid);
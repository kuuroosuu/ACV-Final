% 測試 8K 8-bit wave file reading

clear waveFile
waveFile{1}='Since his wife died, his life has had no meaning_70.wav';
waveFile{2}='lately.wav';

for i=1:length(waveFile)
	fprintf('%s: ', waveFile{i});
	[y1, fs, nbits]=wavread(waveFile{i});
	y1=(y1+1)*2^nbits/2;

	y2=rawRead(waveFile{i});
	y2(1:44)=[];	% 砍掉檔頭部分
	
%	if (nbits==16)
%		y2=y2(1:2:end)+y2(2:2:end)*2^8;
%	end

	fprintf('difference=%g\n', sum(abs(y1-y2)));
end
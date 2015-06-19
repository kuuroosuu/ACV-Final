function resultpitch=pitchLabel(waveFile)
% pitchLabel: Label the pitch of a given wave file

% Roger Jang 20040911

if nargin<1, selfdemo; return; end

switch lower(waveFile(end-2:end))
	case 'wav'
		[y, fs, nbits]=wavread(waveFile);
		wave=y*2^nbits/2;
	case 'raw'
		y=rawread(waveFile);
		wave=y-128;
		fs=8000;
		nbits=8;
	otherwise
		error('Unknown wave file!');
end
PP=setPitchPrm(fs, nbits);
PP.frameSize=fs/8000*256;
PP.overlap=0;
PP.maxShift=PP.frameSize;
PP.waveFile=waveFile;
PP.targetPitchFile=[waveFile(1:end-3), 'Pnsdf'];
%if exist(targetPitchFile)
%	PP.targetPitchFile=targetPitchFile;
%end
resultpitch=wave2pitch(wave, fs, nbits, 1, PP);
pFile=strrep(PP.targetPitchFile, '''', '''''');		% 將「'」帶換成「''」
fp = fopen(PP.targetPitchFile,'w');
fprintf(fp,'%g\r\n',resultpitch);
fclose(fp);
fprintf('Save pitch vector to "%s"\n', PP.targetPitchFile);

% ====== Self demo
function selfdemo
waveFile='倫敦鐵橋垮下來_不詳_0.wav';
feval(mfilename, waveFile);
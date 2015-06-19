function waveEnhancement(wavFile1, wavFile2, plotOpt)
% waveEnhancement: Wave enhancement, with zero-mean
%	Usage: waveEnhacement(wavFile1, wavFile2)
%
%	For example:
%		waveFile1='主人下馬客在船.wav';
%		waveFile2='test.wav';
%		waveEnhancement(waveFile1, waveFile2, fs2, nbits2)
%		dos(['start ', waveFile2]);

%	Roger Jang, 20090206

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

[y, fs, nbits]=wavReadInt(wavFile1);
y=y-mean(y);			% Mean-justification
%y2=SSScalart96(y, fs);		% Spectral substraction, which lead to worse performance for VC
%y2=specSubtract(y, fs);
y2=specSubtractMex(y, fs);
%maxY=max(abs(y)); y=y/maxY;	% Volume normalization

if strcmp(wavFile1, wavFile2)	% 若檔案一樣，稍等片刻，比較穩定
	pause(0.1);
end

%fprintf('Saved %s\n', wavFile2);
wavwrite(y2/(2^nbits/2), fs, nbits, wavFile2);

if plotOpt
	subplot(2,1,1); plot((1:length(y))/fs, y);
	subplot(2,1,2); plot((1:length(y2))/fs, y2);
	drawnow;
	sound(y, fs);
	sound(y2, fs);
end

% ====== Selfdemo
function selfdemo
wavFile1='noisy.wav';
wavFile2='enhanced.wav';
plotOpt=1;
feval(mfilename, wavFile1, wavFile2, plotOpt)
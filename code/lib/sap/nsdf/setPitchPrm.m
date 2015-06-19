function PP=setPitchPrm(fs, nbits)
% setParamPrm: 設定 pitch tracking 的參數

if nargin<1, fs=8000; end
if nargin<2, nbits=8; end

% ====== 固定參數，請勿修改
PP.fs=fs;				% 取樣頻率
PP.nbits=nbits;				% 位元解析度

% 音框長度及相鄰音框重疊點數
switch(PP.fs)
	case 8000	% 配合 Micron 使用
		PP.frameSize=256;
		PP.overlap=0;
	case 16000	% frameSize 是 HTK 的兩倍，搭配語音辨識使用
		PP.frameSize=640;
		PP.overlap=320;
	otherwise	% Take care of fs=11025, 22050, 44100, etc
		PP.frameSize=round(512*PP.fs/11025);
		PP.overlap=round(340*PP.fs/11025);
end

% ====== 音高追蹤的可變參數，可以修改，也可以加入自己定義的參數
PP.fs=fs;
PP.nbits=nbits;
PP.maxShift=PP.frameSize/2;				% ACF/AMDF 的平移量
PP.minPitch = 30;					% 最低音高
PP.maxPitch = 84;					% 最高音高
PP.absVolTh = 2*PP.frameSize*(2^(PP.nbits-8));		% 靜音的絕對門檻值
PP.checkMultipleFreq=1;					% 是否進行倍頻轉換
PP.maxAmdfLocalMinCount=20*PP.frameSize/256;		% 若超過此個數，視為氣音（此數需要配合 fs 及 frameSize 調整）
%PP.frame2pitchFcn='frame2pitch';			% 由音框計算音高的函數
%PP.frame2pitchFcn='frame2pitchByAcf';
PP.frame2pitchFcn='myFrame2pitch';

PP.maxFreq=pitch2freq(PP.maxPitch);
PP.minFreq=pitch2freq(PP.minPitch);

% ====== 以下是旋律辨識的參數
CP.useRest = 1;				% 是否使用休止符
PP.duration = 5;			% 輸入音訊長度
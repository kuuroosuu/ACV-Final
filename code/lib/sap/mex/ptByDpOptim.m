% Plot the influence of ptOpt.indexDiffWeight on pitch tracking

addMyPath;
clear all; more off; pack

waveFile='d:\users\jang\matlab\toolbox\audioProcessing\山月隨人歸.wav';		% 16 kHz, 16 bits
waveFile='00001.wav';
waveFile='\users\jang\app\asr\source\data\errorFile\ptichError\觀眾掀起哈雷熱.wav';
waveFile='\users\jang\matlab\toolbox\sap\twinkle_twinkle_little_star.wav';
pvFile=[waveFile(1:end-4), '.pv'];

pv=asciiRead(pvFile);
index=find(pv);
pv(pv==0)=nan;
pv=pv';

[y, fs, nbits]=wavReadInt(waveFile);

if fs>16000	% 若 fs 過高，先進行 resample
	fs2=16000;
	y2=resample(y, fs2, fs);
	fs=fs2;
	y=y2;
end
		
pfType=0;	% 0:AMDF, 1:ACF
ptOpt=ptOptSet(fs, nbits, pfType);
ptOpt.pfWeight=1;
weights=linspace(0, 2400, 101);			% amdf, 8k, 8bit
for i=1:length(weights)
	fprintf('%d/%d: indexDiffWeight=%f\n', i, length(weights), weights(i));
	ptOpt.indexDiffWeight=weights(i);
	pitch=wave2pitchByDpMex(y, fs, nbits, ptOpt);
	pitchError1(i)=sum(abs(pv(index)-pitch(index)));
end
subplot(2,1,1); plot(weights, pitchError1, 'b.-'); title('pfType=1 (AMDF)');

pfType=2;	% 1:AMDF, 2:ACF
ptOpt=ptOptSet(fs, nbits, pfType);
ptOpt.pfWeight=1;
weights=linspace(0, 1200, 101);		% amdf, 16k, 16bit
for i=1:length(weights)
	fprintf('%d/%d: indexDiffWeight=%f\n', i, length(weights), weights(i));
	ptOpt.indexDiffWeight=weights(i);
	pitch=wave2pitchByDpMex(y, fs, nbits, ptOpt);
	pitchError2(i)=sum(abs(pv(index)-pitch(index)));
end
subplot(2,1,2); plot(weights, pitchError2, 'b.-'); title('pfType=2 (ACF)');
function pitch=freq2pitch(freq)
% FREQ2PITCH 以查表法計算 frequency 到 pitch (in semitone) 的轉換（for 16-bit MCU）
%	Usage: pitch = freq2pitch(freq)
%		freq: Real freq. multiple by 10
%		pitch: Real semitone multiplied by 10
%	This function is only good for 16-bit MCU. For 8-bit MCU, try "freq2pitch8051.m".

%	Roger Jang, 20020120

if nargin==0, selfdemo; return; end

% 下表之計算方式： i=1:60; pitchTable(i)=round(10*semitone2freq(i+32));
pitchTable = [550 583 617 654 693 734 778 824 873 925 980 1038 1100 1165 1235 1308 1386 1468 1556 1648 1746 1850 1960 2077 2200 2331 2469 2616 2772 2937 3111 3296 3492 3700 3920 4153 4400 4662 4939 5233 5544 5873 6223 6593 6985 7400 7840 8306 8800 9323 9878 10465 11087 11747 12445 13185 13969 14800 15680 16612 
];

if freq<pitchTable(1)
	pitch=0;
	return;
end
for i = 2:length(pitchTable),
	if freq<pitchTable(i),
		gap=pitchTable(i)-pitchTable(i-1);
		pitch=10*(32+i-1)+floor((10*(freq-pitchTable(i-1))+floor(gap/2))/gap);	% Linear interpolation
		return;
	end
end
pitch = 0;

function selfdemo
for freq=1:1800
	pitch(freq)=freq2pitch(freq*10);
end
%pitch(pitch==0)=nan;
plot(pitch, '.-');
xlabel('Frequency (Hz)');
ylabel('Pitch (semitone) x 10');
grid on
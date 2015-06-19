clear all

%fprintf('Compiling epdByVolHodMex.cpp...\n');
%mex epdByVolHodMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp

addMyPath
waveFile='../此恨綿綿無絕期+sil.wav';
waveFile='../九微片片飛花璅.wav';
waveFile='C:\dataSet\資策會-派斌\pmp-epd\0001_16k.wav';
waveFile='../csNthu_omniMic.wav';
waveFile='../漢文皇帝有高臺.wav';
waveFile='../漁翁夜傍西巖宿.wav';
waveFile='../星臨萬戶動.wav';
waveFile='../紫駝之峰出翠釜.wav';
%waveFile='../盤賜將軍拜舞歸.wav';
%waveFile='../長安城頭頭白烏.wav';
%waveFile='../日暮聊為梁父吟.wav';
waveFile='../大智若魚.wav';
waveFile='../問君何所之.wav';
waveFile='../魂返關塞黑.wav';
waveFile='../言入黃花川.wav';
waveFile='../魂返關塞黑.wav';
waveFile='../忽然遭世變.wav';
waveFile='../野哭千家聞戰伐.wav';
waveFile='../漾漾汎菱荇.wav';
waveFile='../我姓山田叫山田純子.wav';
waveFile='C:\dataSet\【資策會】日本華語老師會話錄音\lesson01\student\我是愛知淑得大學的學生.wav';
waveFile='../file19.wav';
waveFile='../file56.wav';
[wave, fs, nbits]=wavReadInt(waveFile);
epdPrm=epdPrmSet(fs);
plotOpt = 1;
[ep1, epFrameId1, segment, vol, hod, vh] = epdByVolHod(wave, fs, nbits, epdPrm, plotOpt);
[ep2, epFrameId2] = epdByVolHodMex(wave, fs, nbits, epdPrm);
subplot(3,1,1);
line(ep2(1)*[1 1]/fs, [min(wave), max(wave)], 'color', 'r', 'linewidth', 2);
line(ep2(2)*[1 1]/fs, [min(wave), max(wave)], 'color', 'r', 'linewidth', 2);
if isequal(ep1, ep2)
	fprintf('Same endpoints!\n');
else
	fprintf('Different endpoints!\n');
	fprintf('Dev1: %g - %g = %g\n', ep1(1), ep2(1), ep1(1)-ep2(1));
	fprintf('Dev2: %g - %g = %g\n', ep1(2), ep2(2), ep1(2)-ep2(2));
end
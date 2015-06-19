clear all;
close all;

% ===== 可調參數設定 =====
fs = 16000;
nbits = 16;
frameSize = 640; % overlap固定為0
frameNum4VolTh = 5; % 開始錄音的前幾個frame為靜音門檻
frameNum4epStart = 3; % 語音前方保留幾個frame的靜音 <--- new
frameNum4seg = 18; % 語音長度最少幾個frame(含開頭及結尾靜音！) <--- modified
frameNum4epEnd = 5; % 靜音超過幾個frame開始設終點 <--- modified
duration = 30; % 錄音最大長度(秒)

% ===== 其他參數設定與計算 =====
AIVoice = analoginput('winsound');
chan = addchannel(AIVoice,1);
set(AIVoice,'SampleRate',fs);
set(AIVoice,'SamplesPerTrigger',fs*duration);
start(AIVoice);

% === 取開始錄音的前"frameNum4VolTh"個frame為靜音門檻
data4VolTh = zeros(frameSize, frameNum4VolTh);
for i=1:frameNum4VolTh
    frame = getdata(AIVoice, frameSize);
    data4VolTh(i,:) = frame;
end
volume = frame2volume(data4VolTh, 1);
volTh = mean(volume);

state = 'readyToSetEpStart';
contFrameNum = 0;
epData = [];
header = zeros(1, frameNum4epStart*frameSize); % <--- new
for i=1:duration*fs/frameSize
    frame = getdata(AIVoice, frameSize);
    vol = frame2volume(frame);
    fprintf('state:%s vol:%f, th:%f\n', state, vol, volTh);

    switch state
        case 'readyToSetEpStart'
            header = [header(frameSize+1:end), frame'];
            if(vol>=volTh)
                state = 'setEpStart';
                epData = [epData header];
            end
        case 'setEpStart'
            epData = [epData frame'];
            if(vol<volTh)
                state = 'readyToSetEpEnd';
            end
        case 'readyToSetEpEnd'
            epData = [epData frame'];
            contFrameNum = contFrameNum + 1;
            if(vol>=volTh && contFrameNum<=frameNum4epEnd)
                state = 'setEpStart';
                contFrameNum = 0;
            end
            if(contFrameNum>frameNum4epEnd)
                state = 'setEpEnd';
            end
        case 'setEpEnd'
            if(length(epData)>=frameNum4seg*frameSize)
                fileName = [strrep(int2str(clock),' ','') '.wav'];
                wavwrite(epData,fs,nbits,fileName);
            end
            contFrameNum = 0;
            epData = [];
            header = []; % <--- new
            state = 'readyToSetEpStart';
        otherwise
            fprintf('!!??\n');
    end
end

stop(AIVoice);
delete(AIVoice);
clear AIVoice;

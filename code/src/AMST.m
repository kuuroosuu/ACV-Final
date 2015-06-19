function AMST( fileName )

addpath ../lib/sap
if ~exist('fileName','var'), fileName = '../media/img2.jpg'; end

% Input score image
oriImg = imread(fileName);
size(oriImg)
imshow(oriImg)

% Convert to binary image
if size(size(oriImg),2) == 3
    biImg = rgb2gray(oriImg);
else
    biImg = oriImg;
end
biImg = biImg>128;

% Find histogram to find horizontal line 
imgHistogram = sum(1-biImg');
plot(imgHistogram)

% Find clusters of set of lines 
x = 1:floor(size(oriImg,1)/10);
gauImgHistogram = filter(gaussmf(x,[floor(size(x,2)/2.5) 0]),1,imgHistogram);
plot(gauImgHistogram)
hold on

% Find cluster center
[pks, locs] = findpeaks(gauImgHistogram);
locs = locs(pks>max(pks)/2);
pks = pks(pks>max(pks)/2);
i = 1;
while 1
    if i>length(locs)-1, break; end
    if locs(i+1)-locs(i) < max(size(oriImg,1)/20,50)
        if pks(i)>pks(i+1)
            rmIdx = i+1;
        else
            rmIdx = i;
        end
        pks(rmIdx) = [];
        locs(rmIdx) = [];
    else
        i = i+1;
    end
end
plot(locs, pks, 'go')
hold off

% Find the space between two clusters
scoreInterval = mean(locs(2:end)-locs(1:end-1));

% Find the lines in the clusters
scoreLines = cell(1,length(locs));
for i = 1 : length(locs)
    sBegin = locs(i)-floor(scoreInterval/2); if sBegin<1, sBegin=1; end
    sEnd = locs(i)+floor(scoreInterval/2); if sEnd>size(oriImg,1), sEnd=size(oriImg,1); end
    [linepks, linelocs] = findpeaks(imgHistogram(1, sBegin:sEnd));
    [linepks, linepksIdx] = sort(linepks, 'descend');
    linelocs = linelocs(linepksIdx);
    for j = 1 : length(linepks)
        % Exclude the lines which are too weak
        if linepks(j)<linepks(1)/1.5 & linepks(j)<size(oriImg,2)*0.8, break; end
        flg = 0;
        for k = 1 : length(scoreLines{i})
            if abs(linelocs(j)-scoreLines{i}(k))<5 % Exclude the lines too close
                flg = 1;
                break
            end
        end
        if ~flg, scoreLines{i} = [scoreLines{i} linelocs(j)]; end
    end
    scoreLines{i} = sort(scoreLines{i}, 'ascend')+sBegin-1;
    if length(scoreLines{i})>5, scoreLines{i} = scoreLines{i}(1:5); end
end

% Find the space between two lines
confirmSpace = 0;
for i = 1 : length(scoreLines)
    if length(scoreLines{i})==5
        confirmSpace = confirmSpace + mean(scoreLines{i}(2:5)-scoreLines{i}(1:4));
        break
    end
end
if ~confirmSpace
    disp('No 5 lines score!')
    for i = 1 : length(scoreLines)
        for j = 1 : length(scoreLines{i})-1
            confirmSpace = min(confirmSpace, scoreLines{i}(j+1)-scoreLines{i}(j));
        end
    end
end

% Mend the lines if the number of lines in the clusters are less than 5
for i = 1 : length(scoreLines)
    if length(scoreLines{i})<5
        lineAddVal = [];
        lineAddPos = [];
        for j = 1 : length(scoreLines{i})
            for k = 1 : 4
                [lnMax, lnMaxIdx] = max(imgHistogram(scoreLines{i}(j)+k*confirmSpace-2) : ...
                    imgHistogram(scoreLines{i}(j)+k*confirmSpace+2));
                lineAddVal = [lineAddVal lnMax];
                lineAddPos = [lineAddPos scoreLines{i}(j)+k*confirmSpace-3+lnMaxIdx];
                [lnMax, lnMaxIdx] = max(imgHistogram(scoreLines{i}(j)-k*confirmSpace-2) : ...
                    imgHistogram(scoreLines{i}(j)-k*confirmSpace+2));
                lineAddVal = [lineAddVal lnMax];
                lineAddPos = [lineAddPos scoreLines{i}(j)-k*confirmSpace-3+lnMaxIdx];
            end
        end
        [lineAddVal, lineAddValIdx] = sort(lineAddVal, 'descend');
        lineAddPos = lineAddPos(lineAddValIdx);
        for j = 1 : length(lineAddVal)
            flg = 0;
            for k = 1 : length(scoreLines{i})
                if abs(lineAddPos(j)-scoreLines{i}(k))<confirmSpace*1.01
                    flg = 1;
                    break
                end
            end
            if ~flg, scoreLines{i} = [scoreLines{i} lineAddPos(j)]; end
            if length(scoreLines{i})==5, break; end
        end
    end
end
imshow(oriImg)
hold on
for i = 1 : length(scoreLines)
    lines = scoreLines{i};
    for j = 1 : length(lines)
        line([1 size(oriImg,2)],[lines(j) lines(j)])
    end
end
hold off

% Remove score line
imgAfterRemoveScoreLine = biImg;
for i = 1 : length(scoreLines)
    for j = 1 : 5
        notePreserve = biImg(scoreLines{i}(j)-2:scoreLines{i}(j)+2,:);
        notePreserve = sum(notePreserve);
        notePreserve = notePreserve>=3;
        rep = zeros(5,size(notePreserve,2));
        for k = 1 : 5
            rep(k,:) = notePreserve;
        end
        imgAfterRemoveScoreLine(scoreLines{i}(j)-2:scoreLines{i}(j)+2,:) = rep;
    end
end
imshow(imgAfterRemoveScoreLine)
% pause

% Find the vertical lines
verLineSet = cell(1,length(locs));
for i = 1 : length(locs)
    sBegin = locs(i)-floor(scoreInterval/2); if sBegin<1, sBegin=1; end
    sEnd = locs(i)+floor(scoreInterval/2); if sEnd>size(oriImg,1), sEnd=size(oriImg,1); end
    verLineHis = sum(1-imgAfterRemoveScoreLine(sBegin:sEnd,:));
    verLine = zeros(1,size(imgAfterRemoveScoreLine,2));
    for j = 1 : size(imgAfterRemoveScoreLine,2)
% Find vertical line position in the histogram
        if verLineHis(j)>3.5*confirmSpace, verLine(j) = 1; end
    end
    
% Plotting found position
%     plot(verLineHis)
%     hold on
%     tverLine = find(verLine==1);
%     tverLine
%     for j = 1 : length(tverLine)
%         plot(tverLine(j),verLineHis(tverLine(j)),'go')
%     end
%     hold off
%     pause
    
% Let the value of histogram of non-vertical-line be zero
    tverLineHis = verLineHis;
    tverLineHis(~logical(verLine)) = 0;
%     plot(tverLineHis)
%     pause
    trueVerLine = zeros(1,size(imgAfterRemoveScoreLine,2));
    while 1
% Trimming the vertical lines that are too close to each other
        [mverLineVal, mverLineValIdx] = max(tverLineHis);
        if ~mverLineVal, break; end
        trueVerLine(mverLineValIdx) = 1;
        tverLineHis(mverLineValIdx-5:mverLineValIdx+5) = 0;
    end
    verLine = find(trueVerLine==1);
    
% Plotting the vertical lines after trimming on histogram
%     plot(verLineHis)
%     hold on
%     verLine
%     for j = 1 : length(verLine)
%         plot(verLine(j),verLineHis(verLine(j)),'go')
%     end
%     hold off
%     pause

% Plotting the vertical lines after trimming on image
%     imshow(biImg)
%     hold on
%     for j = 1 : length(verLine)
%         line([verLine(j) verLine(j)],[sBegin sEnd]) 
%     end
%     hold off
%     pause

% Find the real (continuous) part of vertical lines
    verLinePos = struct;
    for j = 1 : length(verLine)
        verLinePos(j).x = verLine(j);
        verLinePos(j).begin = 1;
        verLinePos(j).end = 1;
        localStart = 1;
        localEnd = 1;
        gap = 0;
        lineImg = imgAfterRemoveScoreLine(sBegin:sEnd,verLine(j));
%         shw = uint8(imgAfterRemoveScoreLine(sBegin:sEnd,verLine(j)-2:verLine(j)+2))*255;
        for k = 2 : sEnd-sBegin+1
%             tshw = shw;
%             tshw(k,3) = 128;
%             imshow(tshw)
%             hold on
            if ~lineImg(k)
                localEnd = localEnd+1;
                gap = 0;
                if localEnd-localStart > verLinePos(j).end-verLinePos(j).begin
                    verLinePos(j).begin = localStart;
                    verLinePos(j).end = localEnd;
                end
            elseif localStart==localEnd
                localStart = k;
                localEnd = k;
                gap = 0;
            else
                gap = gap+1;
            end
% If gap is too large, then consider it is discontinuous
            if gap>2
                localStart = k;
                localEnd = k;
                gap = 0;
            end
%             localStart
%             localEnd
%             verLinePos(j).begin
%             verLinePos(j).end
%             gap
%             line([3 3],[verLinePos(j).begin verLinePos(j).end])
%             hold off
%             drawnow
        end
        verLinePos(j).begin = verLinePos(j).begin + sBegin-1;
        verLinePos(j).end = verLinePos(j).end + sBegin-1;
    end
    
% Removing the vertical line that is too thick
%     noInclude = [];
%     tmp = length(verLinePos);
%     rng = 2;
%     for j = 1 : length(verLinePos)
%         left = max(1,verLinePos(j).x-rng);
%         right = min(size(biImg,2),verLinePos(j).x+rng);
%         imshow(biImg(verLinePos(j).begin:verLinePos(j).end,left:right))
%         sum(sum(~biImg(verLinePos(j).begin:verLinePos(j).end,left:right)))/((2*rng+1)*(verLinePos(j).end-verLinePos(j).begin+1))
%         pause
%         if sum(sum(~biImg(verLinePos(j).begin:verLinePos(j).end,left:right)))/(11*(verLinePos(j).end-verLinePos(j).begin+1))>0.8
%             noInclude = [noInclude j];
%         end
%     end
%     verLinePos(noInclude) = [];
%     if length(verLinePos)<tmp, disp('trim');end

    verLineSet{i} = verLinePos;
end

% Plotting the final vertical lines
for i = 1 : length(verLineSet)
    for j = 1 : length(verLineSet{i})
        line([verLineSet{i}(j).x verLineSet{i}(j).x],[verLineSet{i}(j).begin verLineSet{i}(j).end])
    end
end

% Finding if two notes are connected
connectedPartSet = cell(1,length(verLineSet));
% Initialization
for i = 1 : length(verLineSet)
    for j = 1 : length(verLineSet{i})
        connectedPartSet{i}(j).lt = 0;
        connectedPartSet{i}(j).rt = 0;
        connectedPartSet{i}(j).lb = 0;
        connectedPartSet{i}(j).rb = 0;
    end
end
for i = 1 : length(verLineSet)
    for j = 1 : length(verLineSet{i})
%         Having connection with right note?
        if j < length(verLineSet{i})
            BW = 1-imgAfterRemoveScoreLine(min(verLineSet{i}(j).begin,verLineSet{i}(j+1).begin)-2:max(verLineSet{i}(j).end,verLineSet{i}(j+1).end)+2, verLineSet{i}(j).x:verLineSet{i}(j+1).x);
            CC = bwconncomp(BW);
            isConnected = 0;
%             One connected component usually has connection
            if CC.NumObjects == 1
               isConnected = 1;
            else
%                 Or has strong enough horizontal line
                horHist = sum(BW,2);
                if sum(horHist>size(BW,2)*0.9)
                    isConnected = 1;
                end
            end
            if isConnected
                BW = BW(:, ceil(size(BW,2)*(1/2-1/5)):floor(size(BW,2)*(1/2+1/5)));
                topPart = sum(BW(1:round(size(BW,1)/2),:));
                downPart = sum(BW)-topPart;
                if topPart>downPart
%                     Top connection
                    connectedPartSet{i}(j).rt = 1;
                    connectedPartSet{i}(j+1).lt = 1;
                else
%                     Button connection
                    connectedPartSet{i}(j).rb = 1;
                    connectedPartSet{i}(j+1).lb = 1;
                end
            end
        end
    end
end


% Finding note region aside vertical line
halfSpace = round(confirmSpace/2);
regionSize = confirmSpace*confirmSpace;
noteRegion = cell(1,length(verLineSet));
for i = 1 : length(verLineSet)
    noteRegion{i} = [];
    for j = 1 : length(verLineSet{i})
%         Region of four corners
%         whole = imgAfterRemoveScoreLine(verLineSet{i}(j).begin-halfSpace:verLineSet{i}(j).end+halfSpace,verLineSet{i}(j).x-confirmSpace:verLineSet{i}(j).x+confirmSpace);
        r1 = imgAfterRemoveScoreLine(verLineSet{i}(j).begin-halfSpace:verLineSet{i}(j).begin+halfSpace,verLineSet{i}(j).x-confirmSpace:verLineSet{i}(j).x);
        r2 = imgAfterRemoveScoreLine(verLineSet{i}(j).begin-halfSpace:verLineSet{i}(j).begin+halfSpace,verLineSet{i}(j).x:verLineSet{i}(j).x+confirmSpace);
        r3 = imgAfterRemoveScoreLine(verLineSet{i}(j).end-halfSpace:verLineSet{i}(j).end+halfSpace,verLineSet{i}(j).x-confirmSpace:verLineSet{i}(j).x);
        r4 = imgAfterRemoveScoreLine(verLineSet{i}(j).end-halfSpace:verLineSet{i}(j).end+halfSpace,verLineSet{i}(j).x:verLineSet{i}(j).x+confirmSpace);
%         subplot(3,3,1);imshow(r1)
%         subplot(3,3,3);imshow(r2)
%         subplot(3,3,5);imshow(whole)
%         subplot(3,3,7);imshow(r3)
%         subplot(3,3,9);imshow(r4)
%         pause
%         Four vertexs of region
        r1p = [verLineSet{i}(j).begin-halfSpace verLineSet{i}(j).begin+halfSpace verLineSet{i}(j).x-confirmSpace verLineSet{i}(j).x];
        r2p = [verLineSet{i}(j).begin-halfSpace verLineSet{i}(j).begin+halfSpace verLineSet{i}(j).x verLineSet{i}(j).x+confirmSpace];
        r3p = [verLineSet{i}(j).end-halfSpace verLineSet{i}(j).end+halfSpace verLineSet{i}(j).x-confirmSpace verLineSet{i}(j).x];
        r4p = [verLineSet{i}(j).end-halfSpace verLineSet{i}(j).end+halfSpace verLineSet{i}(j).x verLineSet{i}(j).x+confirmSpace];
        noteRegCan = [r1p; r2p; r3p; r4p];
        regStrength = [mean(mean(r1)) mean(mean(r2)) mean(mean(r3)) mean(mean(r4))];
        if connectedPartSet{i}(j).lt, regStrength(1) = inf; end
        if connectedPartSet{i}(j).rt, regStrength(2) = inf; end
        if connectedPartSet{i}(j).lb, regStrength(3) = inf; end
        if connectedPartSet{i}(j).rb, regStrength(4) = inf; end
        [minR, minRIdx] = min(regStrength);
        if minR<0.5
            noteRegion{i} = [noteRegion{i}; noteRegCan(minRIdx,:)];
            
%             else is bar
        end
    end
end

% Show the result after note finding
imshow(biImg)
hold on
for i = 1 : length(verLineSet)
    for j = 1 : size(noteRegion{i},1)
        y = (noteRegion{i}(j,1)+noteRegion{i}(j,2))/2;
        x = (noteRegion{i}(j,3)+noteRegion{i}(j,4))/2;
        plot(x,y,'ro')
    end
end
hold off

% Make music sound
pitch = [];
duration = [];
for i = 1 : length(verLineSet)
    baseLine = scoreLines{i}(1);
    for j = 1 : size(noteRegion{i},1)
        y = (noteRegion{i}(j,1)+noteRegion{i}(j,2))/2;
        pitch = [pitch 64+round((-y+baseLine)/(confirmSpace/2))];
        if connectedPartSet{i}(j).lt || connectedPartSet{i}(j).rt || connectedPartSet{i}(j).lb || connectedPartSet{i}(j).rb
            duration = [duration 1];
        else
            duration = [duration 2];
        end
    end
end
note = [pitch; duration]; note=note(:)';
wave = note2wave(note, 1, 16000);
wavwrite(wave, 16000, 16, 'result.wav')

end


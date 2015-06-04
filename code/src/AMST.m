function [ output ] = AMST( fileName )

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
pause

% Find the vertical lines
verLineSet = cell(1,length(locs));
for i = 1 : length(locs)
    sBegin = locs(i)-floor(scoreInterval/2); if sBegin<1, sBegin=1; end
    sEnd = locs(i)+floor(scoreInterval/2); if sEnd>size(oriImg,1), sEnd=size(oriImg,1); end
    verLineHis = sum(1-biImg(sBegin:sEnd,:));
    verLine = zeros(1,size(biImg,2));
    for j = 1 : size(biImg,2)
        if verLineHis(j)>3*confirmSpace, verLine(j) = 1; end
    end
    
%     plot(verLineHis)
%     hold on
%     tverLine = find(verLine==1);
%     tverLine
%     for j = 1 : length(tverLine)
%         plot(tverLine(j),verLineHis(tverLine(j)),'go')
%     end
%     hold off
%     pause
    
    tverLineHis = verLineHis;
    tverLineHis(~logical(verLine)) = 0;
%     plot(tverLineHis)
%     pause
    trueVerLine = zeros(1,size(biImg,2));
    while 1
        [mverLineVal, mverLineValIdx] = max(tverLineHis);
        if ~mverLineVal, break; end
        trueVerLine(mverLineValIdx) = 1;
        tverLineHis(mverLineValIdx-5:mverLineValIdx+5) = 0;
    end
    verLine = find(trueVerLine==1);
    
%     plot(verLineHis)
%     hold on
%     verLine
%     for j = 1 : length(verLine)
%         plot(verLine(j),verLineHis(verLine(j)),'go')
%     end
%     hold off
%     pause
    globalStart = sBegin;
    globalEnd = sBegin;
    globalCnt = 0;
    localStart = sBegin;
    localEnd = sBegin;
    localCnt = 0;
    for j = 1 : length(verLine)
        for k = sBegin : sEnd
            if ~biImg(k,verLine(j))
                
            end
        end
    end
end

output = scoreLines;

end


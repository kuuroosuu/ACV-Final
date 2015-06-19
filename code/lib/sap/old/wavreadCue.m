function [cueChk, msg] = wavreadCue(file, fs)

chk=[99 117 101 32];    %Cue Chunk ID
fid=fopen(file);
temp=fread(fid);
fclose(fid);
cueStr=temp';


curCuePoint = find(cueStr==chk(1));
inx=[];
for i=1:length(curCuePoint)-1
    chkNum = 0;
    for j=2:4
        if cueStr(curCuePoint(i)+j-1)==chk(j),
            chkNum = chkNum + 1;
        end
    end
    if chkNum == 3,
        inx = [inx curCuePoint(i)];
    end
end

tapNum = cueStr(inx+8)-11;    %get the total number of cue chunk

for i=1:tapNum
    tmp=[];
    for j=4:-1:1
       if length(num2str(cueStr(inx+16+24*(i-1)+j-1)))==1
           tmp = [tmp '0' dec2hex(cueStr(inx+16+24*(i-1)+j-1))];
       else
           tmp = [tmp dec2hex(cueStr(inx+16+24*(i-1)+j-1))];
       end
   end
   tapDuration(i) = hex2dec(tmp);
end

if isempty(inx),
    cueChk=[];
    msg='Not a cue chunk of wave';
else
    a=sort(tapDuration);
    cueChk = a/fs;
    msg='correct';
end
function [pitch, index1, selectIndex, clarity] = myNsdf2pitch(nsdf, fs, plotOpt) 

%for i=1:length(nsdf)
%	if nsdf(i)<=0 
%		nsdf(i)=0;
%	end
%end	
maxList=[];
maxValue=[];
tempIndex = 0;
tempMax = 0;
for i=2:length(nsdf)-1
	if nsdf(i) > 0	
		if nsdf(i)>=nsdf(i+1) & nsdf(i)>=nsdf(i-1) & nsdf(i)>tempMax
			tempMax = nsdf(i);
			tempIndex = i;
		end
		if nsdf(i+1) <= 0 & tempIndex ~= 0
			maxList = [maxList tempIndex];
			maxValue = [maxValue tempMax];
			tempMax = 0;
			tempIndex = 0;
		end
	end
end	

[maxKeyValue,index] = max(maxValue);
maxKeyLag = maxList(index);

selectIndex = 1;
constK = 0.8;

for i=1:length(maxValue)
	if maxValue(i) >= constK*maxKeyValue 
		selectIndex = maxList(i);
		break;
	end
end	
if selectIndex == 0
	maxList(i)
	i
	maxList
end
index1 = 1;
if selectIndex ~= 1
	freq=fs/(selectIndex-index1);
	pitch = freq2pitch(freq);
	clarity = nsdf(selectIndex);
else
	pitch = 0;
	clarity = 0;
end	


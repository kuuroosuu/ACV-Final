function out=pitchSmooth(pitch)

out=pitch;
% 0abc ===> 000c if abs(a-b)<1 && abs(b-c)>3
% abc0 ===> a000 if abs(b-c)<1 && abs(a-b)>3
th=3;
for i=2:length(pitch)-2		% Hot spot is 2 of 0abc
	if pitch(i-1)==0 && pitch(i)~=0 && pitch(i+1)~=0 && pitch(i+2)~=0
		if abs(pitch(i)-pitch(i+1))<1 && abs(pitch(i+1)-pitch(i+2))>th
			out(i)=0; out(i+1)=0;
		end
	end
	if pitch(i-1)~=0 && pitch(i)~=0 && pitch(i+1)~=0 && pitch(i+2)==0
		if abs(pitch(i)-pitch(i+1))<1 && abs(pitch(i)-pitch(i-1))>th
			out(i)=0; out(i+1)=0;
		end
	end
end
% 處理頭尾
if prod(pitch(1:3))~=0 && abs(pitch(1)-pitch(2))<1 && abs(pitch(2)-pitch(3))>th
	out(1)=0; out(2)=0;
end
if prod(pitch(end-2:end))~=0 && abs(pitch(end)-pitch(end-1))<1 && abs(pitch(end-1)-pitch(end-2))>th
	out(end)=0; out(end-1)=0;
end

pitch=out;
% 0ab ===> 00b if abs(a-b)>10
% ab0 ===> a00 if abs(a-b)>10
th=2;
for i=2:length(pitch)-1
	if pitch(i-1)==0 && abs(pitch(i)-pitch(i+1))>th, out(i)=0; end
	if pitch(i+1)==0 && abs(pitch(i)-pitch(i-1))>th, out(i)=0; end
end

% Repeat the previous one
pitch=out;
th=2;
for i=2:length(pitch)-1
	if pitch(i-1)==0 && abs(pitch(i)-pitch(i+1))>th, out(i)=0; end
	if pitch(i+1)==0 && abs(pitch(i)-pitch(i-1))>th, out(i)=0; end
end

% Repeat the previous one
pitch=out;
th=2;
for i=2:length(pitch)-1
	if pitch(i-1)==0 && abs(pitch(i)-pitch(i+1))>th, out(i)=0; end
	if pitch(i+1)==0 && abs(pitch(i)-pitch(i-1))>th, out(i)=0; end
end

pitch=out;
% abc ===> axc where x=(a+c)/2, if b-a>10 & b-c>10
th=3;
for i=2:length(pitch)-1
	if pitch(i-1)>0 && pitch(i+1)>0
		if pitch(i)-pitch(i-1)>th && pitch(i)-pitch(i+1)>th, out(i)=(pitch(i-1)+pitch(i+1))/2; end
		if pitch(i)-pitch(i-1)<-th && pitch(i)-pitch(i+1)<-th, out(i)=(pitch(i-1)+pitch(i+1))/2; end
	end
end

pitch=out;
% 處理頭尾
if pitch(1)~=0 && pitch(2)==0 && pitch(3)==0
	out(1)=0;
end
if pitch(end-2)==0 && pitch(end-1)==0 && pitch(end)~=0
	out(end)=0;
end

%out(out==0)=nan;
%out=medianFilter(out, 7);
%out(isnan(out))=0;
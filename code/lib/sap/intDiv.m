function out=intDiv(x, y)
% 模擬 C 的整數除法，並經由四捨五入來保留精確度

if x==0,
	out=0;
elseif x>0
	out=fix((x+fix(y/2))/y);
else
	out=fix((x-fix(y/2))/y);
end
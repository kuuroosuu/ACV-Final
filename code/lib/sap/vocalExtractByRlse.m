function extracted=vocalExtractByRlse(x1, x2, opt)
%vocalExtractByRlse: Vocal extraction by RLSE

%	Roger Jang, 20121113

if nargin<3
	opt.order=5
	opt.lambda=0.999;
end

% Reset the static variables in rlseOneUpdate.mex
clear rlseOneUpdate

n=length(x1);
extracted=zeros(n,1);
output=zeros(n,1);
for i=opt.order:n
%	fprintf('%d/%d\n', i, rowNum);
	output(i)=rlseOneUpdate(x1(i-opt.order+1:i), x2(i), opt.lambda);
	extracted(i)=x2(i)-output(i);
end

% Reset the static variables in rlseOneUpdate.mex
clear rlseOneUpdate

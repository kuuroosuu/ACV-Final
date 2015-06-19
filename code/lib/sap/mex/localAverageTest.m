% localAverage.m 和 localAverageMex.dll 的同步測試

mex localAverageMex.cpp

n=10000;
for i=1:n
	fprintf('%d/%d\n', i, n);
	frame=round(100*randn(1, 3));
	frame1=localAverage(frame);
	frame2=localAverageMex(frame);
	if sum(abs(frame1-frame2))~=0
		fprintf('Error!');
		return;
	end
end
n=10000;
for i=1:n
	fprintf('%d/%d\n', i, n);
	frame=randn(320, 1);
	zcr1=zcRate(frame);
	zcr2=zcRateMex(frame);
	if (zcr1~=zcr2)
		break;
	end
end

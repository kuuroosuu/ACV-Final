function out = medianFilter(in, order, plotOpt)
%medianFilter: Median filter for a given filter
%	Usage: out = medianFilter(in, order)

%	Roger Jang, 20040424

if nargin<1, selfdemo; return; end
if nargin<2, order=3; end
if nargin<3, plotOpt=0; end

[m, n]=size(in);

out=in;
side=fix(order/2);
if (m==1)|(n==1)	% input is a vector
	for i=(side+1):length(in)-side
		out(i)=median(in((i-side):(i+side)));
	end
else			% input is a matrix
	for i=1:n
		out(:,i)=medianFilter(out(:,i), order);
	end
end

if plotOpt
	plot(1:length(in), in, '-o', 1:length(out), out, '-o');
	legend('Original', ['Order=', num2str(order)]); 
end

% ====== self demo
function selfdemo
n=20;
in=round(10*rand(1, n));
subplot(2,1,1); out1=feval(mfilename, in, 3, 1);
subplot(2,1,2); out2=feval(mfilename, in, 5, 1);
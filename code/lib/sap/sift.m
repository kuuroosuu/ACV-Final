function [yEstimated, excitation, coef]=sift(y, lpcOrder, plotOpt);

if nargin<1, plotOpt=0; end
if nargin<2, lpcOrder=20; end
if nargin<3, plotOpt=1; end

coef = lpc(y, lpcOrder);
yEstimated = filter([0 -coef(2:end)], 1, y);	% Estimated signal
excitation = y - yEstimated;			% Prediction error

if any(isnan(coef))	% Where nan occurs when y is a zero vector?
	coef=zeros(size(coef));
	yEstimated=zeros(size(yEstimated));
	excitation=zeros(size(excitation));
end
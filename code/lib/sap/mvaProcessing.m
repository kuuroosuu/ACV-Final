function out=mvaProcessing(feaMat, opt)
% mvaProcessing: MVA processing of speech features
%
%	Usage:
%		out=mvaProcessing(feaMat)
%		out=mvaProcessing(feaMat, opt)
%
%	Description:
%		out=mvaProcessing(feaMat, opt) returns the MVA processing on the feature matrix feaMat
%
%	Example:
%		feaFile='test.fea';
%		feaObj=feaRead(feaFile);
%		mfcc=feaObj.feature;
%		mfcc2=mvaProcessing(mfcc);
%		subplot(2,1,1); mesh(mfcc);  daspect([1 1 2]);   box on; rotate3d on
%		xlabel('Frame index'); ylabel('Dimension'); title('Before MVA processing');
%		subplot(2,1,2); mesh(mfcc2); daspect([1 1 0.2]); box on; rotate3d on
%		xlabel('Frame index'); ylabel('Dimension'); title('After MVA processing');
%
% 	Reference:
%		Chia-Ping Chen and J. Bilmes, "MVA Processing of Speech Features", IEEE Transactions on Audio, Speech, and Language Processing, 15 (1), pp. 257-270, January 2007

if nargin<1, selfdemo; return; end
% ====== Set the default options
if ischar(feaMat) & strcmp(lower(feaMat), lower('defaultOpt'))
	out.m=2;
	return
end
if nargin<2, opt=feval(mfilename, 'defaultOpt'); end

[dim, frameNum]=size(feaMat);
out=feaMat;

% ====== Perform mean subtraction and variance normalization
mu=mean(feaMat, 2);
variance=var(feaMat, [], 2);
for i=1:dim
	out(i,:)=(out(i,:)-mu(i))/sqrt(variance(i));
end

% ====== Perform ARMA
if opt.m>1
	out2=out;
	for i=opt.m+1:frameNum-opt.m
		out2(:,i)=mean(out(:, i-opt.m:i+opt.m), 2);
	end
	out=out2;
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);

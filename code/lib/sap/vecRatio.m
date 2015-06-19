function output=vecRatio(vec)
% vecRatio: Convert a vector into its ratio format
%	Usage: output=vecRatio(vec)
%		The output is equal to (vec-min(vec))/(max(vec)-min(vec))

%	Roger Jang, 20101029

vecMin=min(vec);
vecMax=max(vec);
if vecMin==vecMax
	fprintf('Warning: the range of the input vec is zero!\n');
	output=0.5*ones(size(vec));
	return
end

output=(vec-vecMin)/(vecMax-vecMin);
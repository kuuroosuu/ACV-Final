function odPrm=odPrmSet
% odPrmSet: Set parameters for OD (onset detection)
%	Usage: odPrm=odPrmSet

%	Roger Jang, 20090608

odPrm.frameSize=0.02;			% In seconds
odPrm.frameShift=1/1000;		% Keep the frame rate to 1000
odPrm.volRatio=0.5;
odPrm.maxTappingPerSec=10;

odPrm.useHighPassFilter=0;
odPrm.filterOrder=3;
odPrm.cutOffFreq=100;
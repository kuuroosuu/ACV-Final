function pitchPrm = setPP(VCfeatures)

%SETPP Set Pitch Parameters.
%   The program is only for VC related codes.
%
%   Cheng-Yuan Lin, 2003, January, 10.

%These following codes are mainly for ptrackfcn function.
%至少一個音框裡要兩個週期(以最小frequency = 50計算, fs = 8000點時)
pitchPrm.frameSize = 512;
pitchPrm.overlap = 0;
pitchPrm.maxFreq = 600;
pitchPrm.lowpassOrder = 5;
pitchPrm.clipThred = 0.3;
pitchPrm.maxThred = 0.1;
pitchPrm.ptrakstyle = 0; %1用frame2pitch, 0用ptrackfcn.
pitchPrm.engthred = 0.2;

%These following codes are mainly for tvc function.
if nargin~=0,  
   pitchPrm.frameSize  = VCfeatures.frameSize;
   pitchPrm.overlap    = VCfeatures.overlap;
   pitchPrm.ptrakstyle = 0; %1用frame2pitch, 0用ptrackfcn.
   %小於此engthred的一定為silence, discard it!
   pitchPrm.engthred = VCfeatures.frameSize*power(0.02,2);
   pitchPrm.zcrthred = VCfeatures.frameSize*(1/20);
   pitchPrm.lowpitch = 60;
   pitchPrm.highpitch = pitchPrm.maxFreq;
end;
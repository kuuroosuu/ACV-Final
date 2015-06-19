function PP=setPitchPrm4dp
% setPitchPrm4dp: 設定 pitch tracking using DP 的參數

%	Roger Jang, 20040524

% ====== Pitch tracking parameters (Used only in wave2pitch.m and midi2mid.m)
PP.ppcNum=5;			% No. of pitch period candidates
PP.ppcIndexDiffWeight=100;	% Weight for PPC index difference between neighboring frames
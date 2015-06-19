function prm=reverbMoorerPrmSet

prm=[];
alpha=5;
% ====== tap delay line
tdl.delay=alpha*[190, 759, 44, 190, 9, 123, 706, 119, 384, 66, 35, 75, 419, 4, 79, 66, 53, 194];
tdl.gain=[.841 .504 .490 .379 .380 .346 .289 .272 .192 .193 .217 .181 .180 .181 .176 .142 .167 .134];
prm=[prm, tdl.delay, tdl.gain];
% ====== lpComb
lpCombPrmDelay=alpha*[1579, 1949, 2113, 2293, 2467, 2647];
lpCombPrmA=0.2*ones(1,6);
lpCombPrmG=[0.1 0.2 0.3 0.2 0.15 0.1];
prm=[prm, lpCombPrmDelay, lpCombPrmA, lpCombPrmG];
% ====== all-pass filter
apfPrm.delay=alpha*307;
apfPrm.bl=0.9;		% apfPrm.fb=-apfPrm.bl; apfPrm.ff=1;
prm=[prm, apfPrm.delay, apfPrm.bl];

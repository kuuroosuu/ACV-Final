function prob=transProb(A, i, k, j)

% 從A(i).candidate(k) 跳到 A(i+1).candidate(j) 的機率

%prob=0; return
%prob=1/abs(A(i).candidate(k).index-A(i+1).candidate(j).index); return

% 根據 pitch 來算 transition probability
pitch1=A(i).candidate(k).pitch;
pitch2=A(i+1).candidate(j).pitch;
width=40;	% prob at x=width is equal to 0.5;
sigma=width/0.8326;
prob=-((pitch1-pitch2)/sigma)^2;	% Log prob. of Gaussian-like PDF

% 根據 amdf 的 local minima 的 index 來算 transition probability
index1=A(i).candidate(k).index;
index2=A(i+1).candidate(j).index;
width=10;	% prob at x=width is equal to 0.5;
sigma=width/0.8326;
prob=-((index1-index2)/sigma)^2;	% Log prob. of Gaussian-like PDF
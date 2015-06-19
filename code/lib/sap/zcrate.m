function out=zcrate(vector)
% zcRate: Zero-crossing rate of a vector

%	Roger Jang, 20050729

%out=sum(abs(diff(vector>0)));			% 此列容易讓靜音的過零率偏高
out=sum(vector(1:end-1).*vector(2:end)<0);	% 位於零點上的資料點，不算過零！

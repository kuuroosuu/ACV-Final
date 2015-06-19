function fpft=fastpartialfouriertransform(x,n,m,shift);

% FPFT fast partial fourier transform
%
% FPFT(X) is equal to FFT(X)
% FPFT(X,N) is equal to FFT(X,N), X is padded with zeroes to length N
% FPFT(X,N,M) is equal to the first M points (rows) of the FFT(X,N)
% FPFT(X,N,M,SHIFT) is equal to the points (rows) SHIFT+1:SHIFT+M of the FFT (X,N)
% N, M and size(X,1) must be a power of 2
%
% For matrices FPFT (like FFT) operates along the 1st dimension of X.

% FPFT was based on the description in
% The Fractional Fourier Transform and Applications
% Daivd H. Bailey and Paul N. Swarztrauber
% October 19, 1995
% SIAM Review, vol. 33 no 3 (Sept. 1991), pg 389-404
% Download: http://crd.lbl.gov/~dhbailey/dhbpapers/

% input argument checking

if (nargin<1)
    error('Not enough input arguments.')
end

if (ndims(x)>2)
    error('X must be vector or 2-D matrix.');
end;

if (nargin<2) % just return the fft
    fpft=fft(x);
    return;
end;

if (nargin<3 || isempty(m) || m==n) % just return the fft
    fpft=fft(x,n);
    return;
end;

% vec=0 array, vec=1 column vector, vec=2 row vector

[lx wx]=size(x);

if wx~=1
    if lx~=1
        vec=0; % 2 dim array
    else
        x=x.';
        lx=wx;
        wx=1;
        vec=2; % row vector
    end
else % column vector
    vec=1; 
end;

if isempty(n)
    n=lx;
else
    if (n<lx),  % truncate x
        x=x(1:n,:);
        lx=n;
    end;
end;

if (m>n)
    error('M must be <= N.');
end;

if rem(log2(n),1) || rem(log2(m),1) || rem(log2(lx),1)
    error('N, size(X,1) and M must be powers of 2.');
end

if nargin<4 || isempty(shift)
    shift=0;
else
    shift=shift-round(shift/n)*n; % shift must be -n/2 <= shift <= n/2
end;

if  m>=n/2 || lx>=n/4 % in this case just do a (more efficient) normal FFT
    fpft=fft(x,n);
    
    % this applies the circular shift
    
    range=[min(n+1,n+shift+1):min(n,n+shift+m)   max(1,shift+1):min(shift+m,n) max(1,shift+m-n):max(0,shift+m-n) ];
    
    switch vec,
        case 0
            fpft=fpft(range,:);
        case 1
            fpft=fpft(range);
        case 2;
            fpft=fpft(range).';
    end;
    
    return;
end;

% now the real magic begins
% from here use the fast fractional fourier transform

persistent expf lookupn countn

maxn=10;  % maximum number of cachable n's, increase this number if needed

if (isempty(countn)) % initialize persistent variables
    lookupn=repmat([NaN 0],maxn,1);   % cache hit counter
    expf=zeros(maxn,n+1);
    countn=0;
end;

indn=find(lookupn(:,1)==n);

if (isempty(indn))
    if (countn>=maxn), % no more cache-able
        [mn,indn]=min(lookupn(:,2));  % delete least used
    else
        countn=countn+1;
        indn=countn;
    end;     
    lookupn(indn,:)=[n 0];
    expf(indn,1:n+1)=exp(-i*pi/n*(0:n).^2);
else
    lookupn(indn,2)=lookupn(indn,2)+1; % increase usage count
end;

indwx=ones(1,wx);
indnm=indn(indwx);

ffty=fft(x.*expf(indnm,1:lx).',2*lx);

if (m>lx)
    fpft=zeros(wx,m);
    
    for s=shift:lx:shift+m-1
        indz=abs([s:lx-1+s -lx+s:-1+s])+1;
        ftz=fft(1./expf(indn,indz).');
        cv=ifft(ftz(:,indwx).*ffty).';
        fpft(:,s-shift+1:s-shift+lx)=cv(:,1:lx).*expf(indnm,indz(1:lx));
    end;
else
    indz=abs([shift:lx-1+shift -lx+shift:-1+shift])+1;
    ftz=fft(1./expf(indn,indz).');
    cv=ifft(ftz(:,indwx).*ffty).';
    fpft=cv(:,1:m).*expf(indnm,indz(1:m));
end;

if (vec~=2)
    fpft=fpft.';
end;

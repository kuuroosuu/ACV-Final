function [y,ff,tt] = specgram(x,n,sr,w,ov,mode,padding)

if (size(x,1) > size(x,2))
  x = x';
end
if (size(x,1)~=1)
    x = x(1,:);
end
% 璝ㄓ方癟痻皚计秖︽计秖玥锣竚

s = length(x);

if nargin < 2
  n = 256;
end
% 璝块セㄧΑま计计秖2(倒x)玥砞﹚ n=256
if nargin < 3
  sr = 1;
end
% 璝块セㄧΑま计计秖3(倒x籔n)玥砞﹚ sr=1
if nargin < 4
  w = n; % window 
end
% 璝块セㄧΑま计计秖4(倒x籔n籔sr)玥砞﹚ w=n
if nargin < 5
  ov = w/2; % overlap
end
% 璝块セㄧΑま计计秖5玥砞﹚ ov=w/2

h = w - ov;

halflen = w/2; % window 
halff = n/2;   % window い翴
acthalflen = min(halff, halflen);

halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen)); % ノcosㄧ计window
win = zeros(1, n);
% ノwindowЧ俱window
win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);

switch mode
    case 0
        %[d ff tt] = spectrogram(x,win,ov,w,sr);
        % 箇皌竚块皚

        ncols = 1+fix((s-n)/h);
        d = zeros((1+n/2), ncols);

        c = 1;

        for b = 0:h:(s-n)
          u = win.*x((b+1):(b+n)); % р块癟–nsamplewindow
          t = fft(u); % 暗FFT
          d(:,c) = t(1:(1+n/2))';
          c = c+1;
        end;

        if (padding)
            tt = [0:h:(s-n)]/sr;
        else
            tt = [n/2:h:(s-n/2)]/sr;
        end
        ff = [0:(n/2)]*sr/n;

    case 1
        %[d ff tt] = spectrogram(x,win,ov,60*w,sr);

        ncols = 1+fix((s-n)/h);
        % Change to be 4k by Frank 2011/5/2
        fb = 32*8/(sr/n);
        %fb = 32*4/(sr/n);
        d = zeros(fb, ncols);

        c = 1;

        for b = 0:h:(s-n)
          u = win.*x((b+1):(b+n)); % р块癟–nsamplewindow
          t = fpft(u,32*n,fb); % Partial FFT : It is even faster if the signal is sparse.
          d(:,c) = t';
          c = c+1;
        end;

        if (padding)
            tt = [0:h:(s-n)]/sr;
        else
            tt = [n/2:h:(s-n/2)]/sr;
        end
        ff = [0:fb-1]*sr/(32*n);

end


if nargout < 1
  imagesc(tt,ff,20*log10(abs(d)));
  axis xy
  xlabel('Time / s');
  ylabel('Frequency / Hz');
else
%if nargout < 1
%   figure;
%   subplot(311);
%   db=20*log10(abs(d));
%   imagesc(tt,ff*1000,db); axis normal;
%   colorbar;
%   axis xy
%   xlabel('Time / s');
%   ylabel('Frequency / Hz');
%   cmap=colormap;
%else
  y = d;
end

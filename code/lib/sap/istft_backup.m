function x=istft(fftMat, hopSize)
% X = istft(fftMat, F, W, H)                   Inverse short-time Fourier transform.
%	Performs overlap-add resynthesis from the short-time Fourier transform 
%	data in D.  Each column of D is taken as the result of an F-point 
%	fft; each successive frame was offset by H points. Data is 
%	hann-windowed at W pts.  
%       W = 0 gives a rectangular window; W as a vector uses that as window.

[dim, frameNum]=size(fftMat);
frameSize=(dim-1)*2;
xlen=frameSize+(frameNum-1)*hopSize;
x=zeros(1, xlen);

win=sin(pi*[0:frameSize-1]/(frameSize-1));	% sine window
for i=1:frameNum
	ft = fftMat(:, i).';
	ft=[ft, conj(ft([((frameSize/2)):-1:2]))];
	px=real(ifft(ft));
	index1=(i-1)*hopSize+1;
	index2=index1+frameSize-1;
	x(index1:index2)=x(index1:index2)+px.*win;
end

% This script is called by goDesignClassifier.m and goTest.m

% Plot 2D Gaussian PDF
figure
pointNum=100;
surfObj=qcSurface(DS2, pointNum, qcParam);		% Compute the surface for each Gaussian PDF
classNum=length(surfObj.class);
% Plot the surface
for i=1:classNum
	subplot(2,2,i);
	mesh(surfObj.xx, surfObj.yy, surfObj.class(i).surface);
	title(['PDF of class ', num2str(i)]);
	axis([-inf inf -inf inf -inf inf]);
end
% Plot the contours
for i=1:classNum
	subplot(2,2,2+i);
	contourf(surfObj.xx, surfObj.yy, surfObj.class(i).surface, 30);
	shading flat; colorbar
end

figure
decisionBoundaryPlot(surfObj);		% Plot the decision boundary
title('Decision boundaries using quadratic classifiers');
qcTrain(DS2, priors, 1);		% Overlay the data points
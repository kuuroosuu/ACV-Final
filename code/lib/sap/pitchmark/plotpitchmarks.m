function plotpitchmarks(y, fs, pitchmarks, SFS_pitchcurve, AXISsetting1, AXISsetting2)

Marks_pitchcurve = (fs./diff(pitchmarks))'; %row-based vector.

if nargin<5,
   AXISsetting1 = 'tight';   
   MAXpitch = max([Marks_pitchcurve; SFS_pitchcurve])+20;  %Unit : Hz
   MINpitch = min([Marks_pitchcurve; SFS_pitchcurve])-20;  %Unit : Hz
   Curvelength = max([length(Marks_pitchcurve) length(SFS_pitchcurve)]);
   AXISsetting2 = [1 Curvelength MINpitch MAXpitch];
end;

figure('position', [-3 43 1280 646]);
% ==================================================================== %
subplot(2,1,1);
plot((1:length(y)), y, '.-');
for i=1:length(pitchmarks),
   line(pitchmarks(i), y(pitchmarks(i)), 'marker', 'o', 'color', 'r', 'linewidth', 1.6);
end
ylabel('Amplitude'); title('Wave form'); legend('Wave','Pitch marks');
axis(AXISsetting1);

% ==================================================================== %
subplot(2,1,2);
if length(SFS_pitchcurve) > length(Marks_pitchcurve),
   SFS_pitchcurve = SFS_pitchcurve(1: length(Marks_pitchcurve));
else
   SFS_pitchcurve = scalevector(SFS_pitchcurve, length(Marks_pitchcurve))'; %row-based vector
end;
A = plot(SFS_pitchcurve, 'ro-'); set(A, 'linewidth', 1.6);
hold on; A = plot(Marks_pitchcurve, 'bx-'); set(A, 'linewidth', 1.6);
RMSE = sqrt(mean(abs(Marks_pitchcurve - SFS_pitchcurve).^2));
title(['RMSE of SFS-pitch and Marks-pitch is: ' num2str(RMSE) ' (Hz).']);
ylabel('Frequency (Hz)'); xlabel('Frame Index (10 ms)');
legend('Pitch curve estimated by using Speech Filing System software (SFS-pitch)', 'Pitch curve estimated directly by using distance of pitch marks (Marks-pitch)', 2);
axis(AXISsetting2);
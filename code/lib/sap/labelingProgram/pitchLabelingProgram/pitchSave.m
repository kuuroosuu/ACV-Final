userDataGet;
pFile=strrep(ptOpt.targetPitchFile, '''', '''''');		% ±N¡u'¡v±a´«¦¨¡u''¡v
%eval(['save ''', pFile, ''' targetPitch -ascii']);
saveCommand=sprintf('asciiWrite(targetPitch, ''%s'')', pFile);
eval(saveCommand);
fprintf('Save pitch vector to "%s"\n', ptOpt.targetPitchFile);
%userDataGet;
pFile=strrep(PP.targetPitchFile, '''', '''''');		% ±N¡u'¡v±a´«¦¨¡u''¡v
eval(['save ''', pFile, ''' targetPitch -ascii']);
fprintf('Save pitch vector to "%s"\n', PP.targetPitchFile);
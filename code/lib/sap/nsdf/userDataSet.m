userDataList;
for i=1:length(userDataVariables)
	if exist(userDataVariables{i})==1
		eval(['userData.', userDataVariables{i}, '=', userDataVariables{i}, ';']);
	end
end
set(0, 'userdata', userData);
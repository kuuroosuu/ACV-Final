userData=get(0, 'userdata');
userDataList;
for i=1:length(userDataVariables)
	if isfield(userData, userDataVariables{i})
		eval([userDataVariables{i}, '=userData.', userDataVariables{i}, ';']);
	end
end
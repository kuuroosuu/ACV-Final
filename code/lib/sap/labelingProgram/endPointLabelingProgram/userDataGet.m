if exist('asdfghjkl'), error('asdfghjkl already exists!'); end
userData=get(gcf, 'userdata');
userDataList;
for asdfghjkl=1:length(userDataVariables)
	if isfield(userData, userDataVariables{asdfghjkl})
		eval([userDataVariables{asdfghjkl}, '=userData.', userDataVariables{asdfghjkl}, ';']);
	end
end
clear asdfghjkl
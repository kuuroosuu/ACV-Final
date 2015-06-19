if exist('asdfghjkl'), error('asdfghjkl already exists!'); end
userDataList;
for asdfghjkl=1:length(userDataVariables)
	if exist(userDataVariables{asdfghjkl})==1
		eval(['userData.', userDataVariables{asdfghjkl}, '=', userDataVariables{asdfghjkl}, ';']);
	end
end
set(gcf, 'userdata', userData);
clear asdfghjkl
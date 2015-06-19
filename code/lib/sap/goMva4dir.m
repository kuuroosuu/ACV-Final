% Run MVA for each .fea file
% Note that the original fea file will be modified!
% Roger Jang, 20111007

dirName='testFeaDir';
feaData=recursiveFileList(dirName, 'fea');
fileNum=length(feaData);

for i=1:fileNum
	fprintf('%d/%d: feaFile=%s\n', i, fileNum, feaData(i).path);
	inputFeaFile=feaData(i).path;
	outputFeaFile=inputFeaFile;
	feaObj=feaRead(inputFeaFile);
	feaObj.feature=mvaProcessing(feaObj.feature, struct('m', 2));
	feaWrite(feaObj, outputFeaFile);
end

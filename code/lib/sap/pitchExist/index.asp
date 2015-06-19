<%@language=JScript%>
<% title="Listing of MATLAB Programs" %>
<!--#include virtual="/jang/include/editfile.inc"-->
<%
function getMFileH1(fileName){
	fso = new ActiveXObject("Scripting.FileSystemObject");
	fid=fso.OpenTextFile(fileName, 1);
	lines=fid.ReadAll();
	fid.Close();
	
	pattern=/^%\s*\w+\s+(.*)/i;
	pattern=/^%\s*(.*)/i;
	lines=lines.split("\n");
	for (i=0; i<lines.length; i++){
		line=lines[i];
		found=line.match(pattern);
		if (found!=null)
			return(RegExp.$1);
	}
	return("");
}
%>
<html>
<head>
	<title><%=title%></title>
	<meta HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=big5">
	<meta HTTP-EQUIV="Expires" CONTENT="0">
</head>

<body background="/jang/graphics/background/yellow.gif">
<h2 align=center><%=title%></h2>
<hr>

<ul>
<p><li>Steps for running the program:
	<ol>
	<li>Change suvParamSet.m, if necessary, to change the parameters for SUV detection.
	<li>goCollectData.m
		<ul>
		<li>Collect the training/test data.
		<li>Store the data at DS.mat.
		</ul>
	<li>goSelectInput.m
		<ul>
		<li>Load DS.mat for input selection using LOO of KNNR.
		<li>Save the indices of the selected inputs at bestInputIndex.mat.
		</ul>
	<li>goPlotData2d.m: Plot the 2d scatter data of the selected inputs.
	<li>goPlotDataAll.m: Plot all combinations of the inputs of the data.
	<li>goQcDesign.m
		<ul>
		<li>Design a quadratic classifier.
		<li>Save the parameters as qcParam.mat.
		<li>Plot the decision boundary if the feature dimension is 2.
		</ul>
	<li>goGmmcDesign.m
		<ul>
		<li>Find the no. of mixtures for each class of GMM.
		<li>Save the parameters as gmmData.mat.
		</ul>
	</ol>
<p><li>List of M files:<p>
<table border=1 align=center>
<tr>
<th>File Name<th>H1 Help<th>File Size
<%
fso = new ActiveXObject("Scripting.FileSystemObject");
fullPath = Server.MapPath(".");
fd = fso.GetFolder(fullPath);
fc = new Enumerator(fd.Files);
for (; !fc.atEnd(); fc.moveNext()){
	fileName=fc.item()+"";
	items=fileName.split(".");
	ext=items[items.length-1];
	if ((ext=="m")|(ext=="M")){
		f = fso.GetFile(fileName);
		Response.write("<tr>");
		Response.write("<td><a href=\"" + f.Name + "\">" + f.Name + "</a></td>");
		Response.write("<td>" + getMFileH1(fileName) + "&nbsp;</td>");
		Response.write("<td>" + f.Size + "</td>");
	}
}
%>
</table>
</ul>

<hr>

<script language="JavaScript">
document.write("Last updated on " + document.lastModified + ".")
</script>

<a href="/jang/sandbox/asp/lib/editfile.asp?FileName=<%=Request.ServerVariables("PATH_INFO")%>"><img align=right border=0 src="/jang/graphics/invisible.gif"></a>
</body>
</html>
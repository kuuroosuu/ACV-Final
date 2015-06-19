<%@language=JScript%>
<!== Roger Jang, 20041120 -->

<%
function getMFileUsage(fileName){
	var i;
	var fso = new ActiveXObject("Scripting.FileSystemObject");
	var fid=fso.OpenTextFile(fileName, 1);
	var contents=fid.ReadAll();
	fid.Close();
	var lines=contents.split("\n");
	var output=new Object();

	// ====== Check if a function
	pattern=/^\s*function\s/i;
	line=lines[0];
	output.isFunction=false;
	if (line.match(pattern) != null)
		output.isFunction=true;

	// ====== Get H1 help
	pattern=/^%\s*.*: (.*)/i;
	pattern=/^%\s*(.*)/i;
	output.h1="";
	for (i=0; i<lines.length; i++)
		if (lines[i].match(pattern) != null){
			output.h1=RegExp.$1;
			break;
		}
	
	// ====== Usage line
	pattern=/^%\s*Usage:\s*(.*)/i;
	output.usage="";
	startLine=-1;
	for (i=0; i<lines.length; i++)
		if (lines[i].match(pattern) != null){
			output.usage=RegExp.$1;
			startLine=i;
			break;
		}

	// ====== Argument documentation
	pattern=/^%\s*(.*)/i;
	output.argumentDoc="";
	if (startLine>=0)
		for (i=startLine+1; i<lines.length; i++)
			if ((lines[i].match(pattern) != null) && (RegExp.$1 != ""))
				output.argumentDoc = output.argumentDoc + "\\n" + RegExp.$1;
			else
				break;

	// ====== Return the final output
	return(output);
}
%>

<%title="Audio Signal Processing Toolbox"%>
<html>
<head>
	<title><%=title%></title>
	<meta HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=big5">
	<meta HTTP-EQUIV="Expires" CONTENT="0">
</head>

<body>
<h2 align=center><%=title%></h2>
<h3 align=center><a target=_blank href="/jang">Roger Jang</a></h3>
<hr>

<!== 簡要說明 -->
<h3><img src="/jang/graphics/dots/redball.gif">Introduction</h3>
<blockquote>
This toolbox provide some utility functions that are commonly used in audio signal processing.
<p>
If you want to download the toolbox, please click <a href="..\audio.rar">here</a>.
<p>
Before using the toolbox, you should add the following toolboxes into your search path:
<a href="../utility">Utility Toolbox</a>. You can add the required path by invoking the script addMyPath.m under MATLAB.
</blockquote>

<!== 列出函數 -->
<%
functionList=new Array();
i=0;
functionList[i++]="freq2pitch.m";
functionList[i++]="pitch2freq.m";
functionList[i++]="pitchPlay.m";
functionList[i++]="frameFlip.m";
functionList[i++]="frame2volume.m";
functionList[i++]="frame2acf.m";
functionList[i++]="frame2smdf.m";
functionList[i++]="frame2acfOverSmdf.m";
functionList[i++]="frame2zcr.m";
functionList[i++]="frame2sampleIndex.m";
functionList[i++]="rawRead.m";
functionList[i++]="wave2formant.m";
functionList[i++]="wave2mfcc12.m";
functionList[i++]="wave2pitch.m";
functionList[i++]="wavReadInt.m";
functionList[i++]="wave2wave.m";
%>
<h3><img src="/jang/graphics/dots/redball.gif">List of Functions</h3>
<table border=1 align=center>
<tr>
<th>Functions<th>H1 Help<th>Usage<th>Description<th>Size
<%
fso = new ActiveXObject("Scripting.FileSystemObject");
fullPath = Server.MapPath(".");
for (i=0; i<functionList.length; i++){
	fileName=fullPath + "/" + functionList[i];
	ext=fso.GetExtensionName(fileName);
	if ((ext=="m") || (ext=="M")){
		f = fso.GetFile(fileName);
		out=getMFileUsage(fileName);
		if (out.isFunction){
			Response.write("\n<tr>\n");
			Response.write("<td><a target=_blank href=\"" + f.Name + "\">" + f.Name + "</a></td>");
			Response.write("<td>" + out.h1 + "&nbsp;</td>");
			doc=out.argumentDoc;
			doc=doc.replace(/"/g, "\\\"");		// 將字串原有的 " 代換成 \"
			doc=doc.replace(/'/g, "\\\'");		// 將字串原有的 ' 代換成 \'
		//	if (doc!="")
		//		Response.write("<td><a href='javascript:alert(\"" + doc + "\")'>" + out.usage + "</a>&nbsp;</td>");
		//	else
		//		Response.write("<td>" + out.usage + "&nbsp;</td>");
			Response.Write("<td>" + out.usage + "&nbsp;</td>");
			Response.write("<td align=right>\n");
			Response.write(doc);
			Response.write("&nbsp</td>");
			Response.write("<td align=right>" + f.Size + "</td>");
		}
	}
}
%>
</table>

<!== 列出底稿 -->
<%
scriptList=new Array();
i=0;
scriptList[i++]="addMyPath.m";
scriptList[i++]="wavReadTest.m";
%>
<h3><img src="/jang/graphics/dots/redball.gif">List of Scripts</h3>
<table border=1 align=center>
<tr>
<th>File Names<th>Help<th>Size
<%
for (i=0; i<scriptList.length; i++){
	fileName=fullPath + "/" + scriptList[i];
	ext=fso.GetExtensionName(fileName);
	if ((ext=="m") || (ext=="M")){
		f = fso.GetFile(fileName);
		out=getMFileUsage(fileName);
		if (!out.isFunction){
			Response.write("\n<tr>\n");
			Response.write("<td><a href=\"" + f.Name + "\">" + f.Name + "</a></td>");
			Response.write("<td>" + out.h1 + "&nbsp;</td>");
			Response.write("<td align=right>" + f.Size + "</td>");
		}
	}
}
%>
</table>

<hr>
<script language="JavaScript">
document.write("Last updated on " + document.lastModified + ".")
</script>

</body>
</html>

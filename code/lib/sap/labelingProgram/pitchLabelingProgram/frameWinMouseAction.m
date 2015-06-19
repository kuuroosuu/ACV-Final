function frameWinMouseAction(action)
% This function is used in "frame2pitch4labeling".

if nargin<1, action='start'; end

switch(action)
	case 'start', % ====== 開啟圖形視窗
		% 設定滑鼠按鈕被按下時的反應動作
		set(gcf, 'WindowButtonDownFcn', [mfilename ' down']);
	case 'down', % ====== 滑鼠按鈕被按下時的反應指令
		% 設定滑鼠移動時的反應指令
		set(gcf, 'WindowButtonMotionFcn', [mfilename, ' move']);
		% 設定滑鼠按鈕被釋放時的反應指令
		set(gcf, 'WindowButtonUpFcn', [mfilename ' up']);
		% 一被按下，即執行 checkPitch move
		feval(mfilename, 'move');
	case 'move', % ====== 滑鼠移動時的反應指令
		% 取得相關資料
		userDataGet;
		acfX=get(acfH, 'xdata');
		acfY=get(acfH, 'ydata');
		% 取得滑鼠座標
		axes(acfAxisH);
		acfXLim=get(acfAxisH, 'xlim');
		currPt=get(gca, 'CurrentPoint');
		mouseX=currPt(1,1);
		mouseY=currPt(1,2);
		% 更新 manualBarH
		if (acfXLim(1)<=mouseX) & (mouseX<=acfXLim(2))
			% 找出距離最近的點
		%	[minDist, pitchIndex]=min(abs(mouseX-acfX));
			localMaxIndex=find(localMaxPos);
			if ~isempty(localMaxIndex)
				[minDist, ind]=min(abs(mouseX-localMaxIndex));
				pitchIndex=localMaxIndex(ind);
			end
		else
			pitchIndex=0;
		end
		set(manualBarH, 'xdata', pitchIndex*[1 1]);
		% 更新 mainWindow 的圖
		if pitchIndex==0,
			targetPitch(frameIndex)=0;
		else
			freq=ptOpt.fs/(pitchIndex-1);
			targetPitch(frameIndex)=freq2pitch(freq);
		end
		temp=targetPitch; temp(temp==0)=nan;
		set(targetPitchH, 'ydata', temp);
		userDataSet;
	case 'up', % ====== 滑鼠按鈕被釋放時的反應指令
		% 清除滑鼠移動時的反應指令
		set(gcf, 'WindowButtonMotionFcn', '');
		% 清除滑鼠按鈕被釋放時的反應指令
		set(gcf, 'WindowButtonUpFcn', '');
end
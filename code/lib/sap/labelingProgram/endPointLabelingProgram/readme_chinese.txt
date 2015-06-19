人工標示端點的程式
Roger Jang, 20050314, 20060304, 20201129

1. 請先下載兩個工具箱：
	Utility Toolbox: http://mirlab.org/jang/matlab/toolbox/utility.rar
	SAP Toolbox:     http://mirlab.org/jang/matlab/toolbox/sap.rar

2. 展開這兩個工具箱，並放到適當目錄。請修改主程式 go.m 以便將這兩個工具箱加入搜尋目錄。

3. 請修改主程式 go.m 的變數 waveDir，將之指到你的錄音檔案所在的目錄，即可執行此程式，開始標示端點。

4. 此程式碼可以開啟每個檔案，使用者可以拖放第一張圖的紅線來標示端點，最後按下 return（或任一個按鍵），即可存檔並切換到下一個檔案。

5. 在標示端點的過程中，程式碼會自動將端點資訊記錄在檔名上，例如 3a.wav 經過標示存檔後，會產生 3a_849_31226.wav 代表端點是 849 及 31226。

6. 若程式有任何執行上的問題，請找張智星。

7. 檔案說明：

go.m
	主程式，用來標示大量語料
 
endPointLabel.m
	用來標示單一語句

userDataGet.m, userDataSet.m, userDataList.m
	用於 endPointLabel.m 的小程式

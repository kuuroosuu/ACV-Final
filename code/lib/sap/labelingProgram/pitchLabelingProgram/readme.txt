This is the pitch labeling program. Please follow the steps to label the pitch.

1. Modify "go.m" to assign the directory of your recordings to the variable "waveDir". 
2. Execute "go" within MATLAB to generate a figure containing the waveform of the first file for pitch labeling.
3. Please maximize the waveform window for easy viewing.
4. You can click on the wrongly identified pitch to show a new figure of ACF and AMDF.
5. You can click on the right pitch point of ACF to correct the pitch in the waveform figure.
6. If there is no pitch at all, you can click outside the axis to make the pitch equal to 0 in the waveform figure.
7. If the volume is small, the pitch will be unstable. Moreover, there is no pitch for unvoiced sounds.
8. If it's hard to identify the pitch point, please make sure your pitch point should generate a continuous curve in the original figure.
9. Computer detected pitch is pitch1, human labeled pitch is pitch2. Before your labeling, pitch2 is equal to pitch1.
10. You can click "Play wave" to play the original recording.
11. You can click "Play pitch" to play the pitch detected by computer.
12. You can click "Play pitch2" to play the pitch labeled by human. The human-labeled pitch contour should be smooth.
13. When you finish labeling a file, press "enter" to save human-labeled pitch and move to the next file.

If you want to label a single file, try "pitchLabel.m".

Pitch labeling requires practice and experiences. If the pitch is not labeled correctly, it will cause problems in performance evaluation for pitch tracking and melody recognition. So make sure you understand what you are doing for every single detail. If you have any questions, please contact with TA as soon as you can.

The following is Chinese translation:

此為人工標示音高的程式，用法如下。

1. 修改 go.m，將 waveDir 指到你之前的錄音檔案所在的目錄。 
2. 執行 go，即可產生一個視窗，對第一個 wave 檔案進行標示 pitch 的工作。
3. 請將此視窗放到最大，以便觀察、點選。
4. 你可以先點選大圖錯誤的音高點，此時會產生小圖，顯示 AMDF 和 ACF。
5. 若在小圖的 ACF 點選最大值的位置（對應到基本週期），大圖的音高點會跟著修正。
6. 若音高不存在，請直接在小圖的圖軸外面點一下，大圖的音高值就會被設定為零。
7. 音量太小時，通常音高就會不穩定，或是根本不存在。遇到氣音時，音高也是不存在。
8. 若不容易看出音高點，請參考左右音框，把握一個原則：音高曲線必須是連續平滑的曲線。
9. pitch1 是電腦辨識的結果，pitch2 是人工標示的結果，尚未進行人工標示前，pitch2 = pitch1。
10. 你可以點選 Play wave 以播放原始錄音檔案。
11. 你可以點選 Play pitch 以播放電腦自動找到的音高曲線。
12. 你可以點選 Play pitch2 以播放修正後的音高曲線，不應該有暴音。
13. 每次標示完一首歌曲，按下 Enter，即自動儲存最後人工標示的 pitch 檔案，並可跳到下一首。

若要測試單一歌曲的音高標示，可以直接執行 pitchLabel.m。

標示音高需要一些經驗和練習，若標示不正確，會造成後續許多測試及校正的問題，所以麻煩各位同學務必要小心進行。
因此若有任何疑問，請盡快直接和助教或是老師聯繫。
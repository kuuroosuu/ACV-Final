% 2D/3D scatter plots for data visualization

addMyPath;
load DS.mat

% ====== Compute the data amount for each class (計算每一個類別的資料量)
figure; dsClassSize(DS, 1);
% ====== Plot class info w.r.t. a single feature (進行類別對特徵的作圖)
figure; dsProjPlot1(DS);
% ====== Project the data into 2D spaces (將資料投影到二度空間)
figure; dsProjPlot2(DS); set(gcf, 'name', 'Original Data');
% ====== Project the data into 3D spaces (將資料投影到三度空間)
figure; dsProjPlot3(DS); set(gcf, 'name', 'Original Data');

% ====== Feature-wise normalization to have zero-mean unity-variance (對每一個特徵進行正規化，使其成為mu=0, sigma=1的高斯密度分佈)
DS.input=inputNormalize(DS.input);
% ====== Project the normalized data into 2D spaces (將資料投影到二度空間)
figure; dsProjPlot2(DS); set(gcf, 'name', 'Normalized Data');
% ====== Project the normalized data into 3D spaces (將資料投影到三度空間)
figure; dsProjPlot3(DS); set(gcf, 'name', 'Normalized Data');
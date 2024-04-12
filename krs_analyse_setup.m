    
%% set-up for the Kernel Representational Similarity Regression 
% =============================================================
clear all; 
close all; clc


% % % (1) set-up / generate the data: 
load '/Volumes/WD2T/Tools/CCP/KRS/Test/HCPtestdata.mat'
Y = Yf; %NEOsumf(:,2);
X = Xf;

%% ========================================================
% generate Y,X data with a specified covariance structure. 
% ========================================================
covdat.vY = 3; %40;
covdat.vX = 10;
covdat.cxx = 0.2;
covdat.cyy = 0.4;
covdat.cxy = 0.2;
covdat.n = 100;
covdat.v = 10;
covdat.defcri = 0.05;
[Y,X] = krs_generate_data(covdat);

% X = [X randn(covdat.n,covdat.vX)]

% %% Inspect the data and filter on distance outliers
% % --------------------------------------------------------
% [fdat,filter,nf] = krs_distfilter(Xf,Yf,'plot');
% Xf = fdat{1};
% Yf = fdat{2};


%% (3) compare the data modalities ('before')
% % --------------------------------------------------------
% krs_plot_compare(Y,X);

% 
% %% (4) develop krs_regression
% % KRS = krs_regression(Xf,Yf);
% 
% % split the P data 
% r = randperm(size(Yf,2));
% ny = 5;
% Yfs = Yf(:,r(1:ny));
% Xfs = Yf(:,r(ny+1:end));
% KRS = krs_regression(Xfs,Yfs);
[stats] = krs_regression_v2(X,Y);


% %% (4) test a model (ADD svregression.. and other methods)
% % --------------------------------------------------------
% tic
% [model] = krs_model(Xf,Yf,'method','ridge');
% toc
% 
% 
% %% (5) check improvement ('after')
% % --------------------------------------------------------
% Xm = model.Yfit;
% krs_plot_compare(Yf,Xm);



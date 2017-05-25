%% Instrument classification

clear; clc; close all;


N = 2^12;
sr = 48000;

filepathsax = 'SMSAdata/sax/';
filepathvio = 'SMSAdata/violin/';
filepathcla = 'SMSAdata/clarinet/';
filepathtru = 'SMSAdata/trumpet/';

filename = ['0000'; '0001'; '0002'; '0003'; '0004'; '0005'; '0006'; '0007'; '0008'; '0009'];

% filter
T = triFilterBank(N, sr);
coef = 13; 
dataAmount = 10; 

instrumentAmount = 4; 
%%
saxdata = [];
saxID = 'sax';
figure; hold on; 
for i = 1:10
[s, ~] = audioread([filepathsax filename(i, :) '.wav']);
mfcc = getMFCC(s, N, T, coef);
saxdata(i,:) = [mfcc];
plot(saxdata(i,:))
end
title('Saxophone MFCC')
%% Clarinet 
clarinetdata = [];
claID = 'cla';
figure; hold on; 
for i = 1:10
[s, ~] = audioread([filepathcla filename(i, :) '.wav']);
mfcc = getMFCC(s, N, T, coef);
clarinetdata(i,:) = [mfcc];
plot(clarinetdata(i,:))
end
title('Clarinet MFCC')
%% Trumpet
trumpetdata = [];
truID = 'tru';
figure; hold on; 
for i = 1:10
[s, ~] = audioread([filepathtru filename(i, :) '.wav']);
mfcc = getMFCC(s, N, T, coef);
trumpetdata(i,:) = [mfcc];
plot(trumpetdata(i,:))
end
title('Trumpet MFCC')
%% Violin
violindata = [];
vioID = 'vio';
figure; hold on; 
for i = 1:10
[s, ~] = audioread([filepathvio filename(i, :) '.wav']);
mfcc = getMFCC(s, N, T, coef);
violindata(i,:) = [mfcc];
plot(violindata(i,:))
end
title('Violin MFCC')

%% Build labels and dataset 
labels = [];

for i = 1:40
   if i < 11
       labels = [labels; saxID];
   elseif i < 21
       labels = [labels; claID];     
   elseif i < 31
       labels = [labels; truID];
   elseif i < 41
       labels = [labels; vioID];
   end
end

% Data X
X = [saxdata; clarinetdata; trumpetdata; violindata];
%% dataset
addpath('Toolboxes/prtools')
warning('off','all'); prwarning(0); prwaitbar off; 
z = prdataset(X, labels);

%% 
%%
figure
scatterd(z)
title('First 2 Cepstrum Coeffients')

%% 

% generate test and training indicies with respect to the class frequencies
part = cvpartition(labels,'HoldOut',0.5); % data divided into 10% test and 90% traning

% Make a training and testing pr data set
pr_X_tr = prdataset(X(part.training,:), labels(part.training));
pr_X_tst = prdataset(X(part.test, :), labels(part.test));

%% Parametric and non-parametric classifier
% Train classifiers
% parametric classifiers
pr_nmsc = nmsc(pr_X_tr); % minimum distance classifier - nmsc


nmsc_w = pr_X_tst * pr_nmsc;
[C_nmsc, r_nmsc] = confmat(nmsc_w);

disp('Confusion matrix for minimum distance classifier')
confmat(nmsc_w)
disp('------------------------------------------------------------------')

%%

pr_ldc = ldc(pr_X_tr);   % linear discriminant analysis - ldc

% Linear Discriminant Analysis
ldc_w = pr_X_tst * pr_ldc;
[C_ldc, r_ldc] = confmat(ldc_w);

disp('Confusion matrix for linear discriminant analysis classifier')
confmat(ldc_w)
disp('------------------------------------------------------------------')


%%

pr_qdc = qdc(pr_X_tr);   % quadratic discriminant analysis - qdc


% Quadratic Discriminant Analysis
qdc_w = pr_X_tst * pr_qdc;
[C_qdc, r_qdc] = confmat(qdc_w);
disp('Confusion matrix for quadratic discriminant analysis classifier')
confmat(qdc_w)
disp('------------------------------------------------------------------')

%% 
% non-parametric classifiers 
[pr_knn,nr_nn_used] = knnc(pr_X_tr);  % k-nearest neighbor classifier - knnc

% K-Nearest Neighbor
knn_w = pr_X_tst * pr_knn;
[C_knn, r_knn] = confmat(knn_w);

disp('Confusion matrix for k-nearest neighbor classifier')
confmat(knn_w)
 % k-NN might be subject to overfitting if the number of neighbors is too low
% The knnc classifier calculates the optimal amount of k

fprintf('Number of k''s in k-nearest neighbors classifier: %1.f\n', nr_nn_used)
disp('------------------------------------------------------------------')

% %% 
% pr_svc = svc(pr_X_tr,'p',2,100);   % support vector machine - svc, polynomial kernel with degree 2, C=100 works good on data
% 
% 
% % Support Vector Machine
% svc_w = pr_X_tst * pr_svc;
% [C_svc, r_svc] = confmat(svc_w);
% 
% disp('Confusion matrix for support vector machine')
% confmat(svc_w)
% disp('------------------------------------------------------------------')
% 
%% 
% Test the classifiers
% apply classifier on test data and calculate confusion matrix

% Table of classification errors 
Data_Class_Accuracy = 100*[1-r_nmsc/part.TestSize; 1-r_ldc/part.TestSize; 1-r_qdc/part.TestSize; 1-r_knn/part.TestSize;];
Classification_Errors = [r_nmsc; r_ldc; r_qdc; r_knn;];
c_names = {'MDC'; 'LDA '; 'QDA '; 'k-NN ';};
disp(table(Classification_Errors, Data_Class_Accuracy, 'RowNames',c_names))


%% Plot classifiers in 2D
close all; 
figure
plot_feature = [1 2];
plot_knn = knnc(z(:,plot_feature));  % KNN
scatterd(z(:,plot_feature)); hold on; 
plotc(plot_knn);
title('Scatter plot of KNN classifier')

figure
plot_ldc = ldc(z(:,plot_feature));   % ldc
scatterd(z(:,plot_feature)); hold on; 
plotc(plot_ldc)
title('Scatter plot of LDA classifier')

figure
plot_qdc = qdc(z(:,plot_feature));   % qdc
scatterd(z(:,plot_feature)); hold on; 
plotc(plot_qdc); 
title('Scatter plot of QDC classifier')

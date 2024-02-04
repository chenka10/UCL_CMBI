%% clearing
clc; clear all; close all;
addpath("BallAndStick\");
addpath("DTI_optim\");
addpath("ZeppelinAndStick\");
addpath("ZeppelinAndStickTortuosity\");
addpath("ZeppelinStickDot\");
addpath("ZeppelinTwoSticks\");
addpath("TensorStickDot\");

%% load ISBI2015 data
selected_voxel = 4;
[Avox,qhat,TE,bvals] = LoadISBI2015Data(selected_voxel);

%% setting names and functions for each model

names = {'DTI',...
    'Ball And Stick',...
    'Zeppelin And Stick',...
    'Zeppelin And Stick (Tortuosity)',...
    'Zeppelin Stick Dot',...
    'Zeppelin Two Sticks',...
    'Tensor Stick Dot'};

fit_funcs = {@FitDTI,...
    @FitBallAndStick,...
    @FitZeppelinAndStick,...
    @FitZeppelinAndStickTortuosity,...
    @FitZeppelinStickDot,...
    @FitZeppelinTwoSticks,...
    @FitTensorStickDot};

compute_funcs = {@ComputeDti,...
    @ComputeBallStick,...
    @ComputeZeppelinStick,...
    @ComputeZeppelinStickTortuosity,...
    @ComputeZeppelinStickDot,...
    @ComputeZeppelinTwoSticks,...
    @ComputeTensorStickDot};

N = numel(names);

params = cell(N,1);
resnorms = zeros(N,1);
models_res = cell(N,1);
AIC_all = zeros(N,1);
BIC_all = zeros(N,1);

%% fit models (and compute AIC/BIC)

K_AIC = numel(Avox);

for i=1:N       
    fit_func = fit_funcs{i};
    [real_params, optim_params, success_rate, min_resnorm] = fit_func(Avox,qhat,bvals);
    params{i} = real_params;
    resnorms(i) = min_resnorm;
    compute_func = compute_funcs{i};
    models_res{i} = compute_func(params{i},bvals,qhat);    

    N_AIC = numel(params{i})+1; % plus one since noise is also estimated
    AIC_all(i) = ComputeAIC(N_AIC,K_AIC,resnorms(i));
    BIC_all(i) = ComputeBIC(N_AIC,K_AIC,resnorms(i));

    WriteLineToCSV({names{i},N_AIC,min_resnorm,AIC_all(i),BIC_all(i),success_rate},['models_results_voxel_' num2str(selected_voxel) '.csv']);

    fprintf('%s; SSD=%d, N=%d, AIC=%d, BIC=%d, success_rate=%d\n',names{i},min_resnorm,N_AIC,AIC_all(i),BIC_all(i),success_rate);
end

%% Visualize fit to data

[~,AIC_sort_indexes] = sort(AIC_all);

figure('Position',[100,-100,700,1600]);
for n=1:N

    i = AIC_sort_indexes(n);

    subplot(8,1,n);
    plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
    hold on;
    plot(models_res{i}, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
    xlabel('Sample num.')
    ylabel('Signal Value')
    legend('Data','Model')
    title(['[' names{i} '] SSD = ' num2str(resnorms(i),4) ', AIC = ' num2str(AIC_all(i),8)])
end

%% Cross validate models

% num folds
F = 6;

cross_corr_results = zeros(N,F);

shuffled_indexes = randperm(numel(Avox));
Avox_shuffled = Avox(shuffled_indexes);
qhat_shuffled = qhat(:,shuffled_indexes);
bvals_shuffled = bvals(shuffled_indexes);

for i=1:N
    fit_func = fit_funcs{i};
    compute_func = compute_funcs{i};
    res = CrossValidateModel(fit_func,compute_func,Avox_shuffled,qhat_shuffled,bvals_shuffled,F);
    cross_corr_results(i,:) = res;
end

%% print cross validation results

for n=1:N
    WriteLineToCSV({cross_corr_results(n,1),cross_corr_results(n,2),cross_corr_results(n,3),cross_corr_results(n,4),cross_corr_results(n,5),cross_corr_results(n,6),mean(cross_corr_results(n,:))},'cross_corr_k_fold_6.csv');
    fprintf('%s: mean cross correlation results = %d\n',names{n},cross_corr_results(n));
end


%%

function results = CrossValidateModel(fit_func,compute_func,Avox,qhat,bvals,fold_num)

results = zeros(fold_num,1);

data_size = numel(Avox);
data_range = 1:data_size;
split_size = floor(data_size/fold_num);
for i=1:fold_num
    train_start = (i-1)*split_size;
    train_end = train_start + split_size;
    test_indexes = logical((data_range>train_start).*(data_range<train_end));
    train_indexes = ~test_indexes;
    [real_params, ~,~,~] = fit_func(Avox(train_indexes),qhat(:,train_indexes),bvals(train_indexes));

    test_res = compute_func(real_params,bvals(test_indexes),qhat(:,test_indexes))';

    results(i) = sum((test_res - Avox(test_indexes)).^2);
    fprintf('cross-validation, iter=%d, SSD=%d\n',i,results(i));
end
end
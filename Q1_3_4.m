%% clearing
clc; clear all; close all;
addpath("BallAndStick\");
addpath("DTI_optim\");
addpath("ZeppelinAndStick\");
addpath("ZeppelinAndStickTortuosity\");
addpath("ZeppelinStickDot\");
addpath("ZeppelinTwoSticks\");
addpath("TensorStickDot\");

%% load data
[Avox,qhat,TE,bvals] = LoadISBI2015Data(2);

%%

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

%% fit models
K_AIC = numel(Avox);

for i=1:N      
    fit_func = fit_funcs{i};
    [real_params, optim_params, success_rate, min_resnorm] = fit_func(Avox,qhat,bvals);
    params{i} = real_params;
    resnorms(i) = min_resnorm;
    compute_func = compute_funcs{i};
    models_res{i} = compute_func(params{i},bvals,qhat);
    fprintf('%s; SSD=%d\n',names{i},min_resnorm);

    N_AIC = numel(params{i});
    AIC_all(i) = ComputeAIC(N_AIC,K_AIC,resnorms(i));
    BIC_all(i) = ComputeBIC(N_AIC,K_AIC,resnorms(i));
end

%% compare given values with model values

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

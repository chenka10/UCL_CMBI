%% clearing
clc; clear all; close all;
addpath("BallAndStick\");
addpath("DTI_optim\");
addpath("ZeppelinAndStick\");
addpath("ZeppelinAndStickTortuosity\");
addpath("ZeppelinStickDot\");

%% load data
[Avox,qhat,TE,bvals] = LoadISBI2015Data(1);

%%
N = 4;

names = {'DTI','Ball And Stick','Zeppelin And Stick', 'Zeppelin Stick Dot'};
fit_funcs = {@FitDTI,@FitBallAndStick,@FitZeppelinAndStick,@FitZeppelinStickDot};
compute_funcs = {@ComputeDti,@ComputeBallStick,@ComputeZeppelinStick,@ComputeZeppelinStickDot};
params = cell(N,1);
resnorms = zeros(N,1);
models_res = cell(N,1);

%% fit models
for i=1:N
    fit_func = fit_funcs{i};
    [real_params, optim_params, success_rate, min_resnorm] = fit_func(Avox,qhat,bvals);
    params{i} = real_params;
    resnorms(i) = min_resnorm;
    compute_func = compute_funcs{i};
    models_res{i} = compute_func(params{i},bvals,qhat);
    fprintf('%s; SSD=%d\n',names{i},min_resnorm);
end

%% compare given values with model values
figure('Position',[100,100,1500,700]);
for i=1:N
    subplot(2,2,i);
    plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
    hold on;
    plot(models_res{i}, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
    xlabel('Sample num.')
    ylabel('Signal Value')
    legend('Data','Model')
    title(names{i})
end

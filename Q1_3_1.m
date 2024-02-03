%% clearing
addpath("BallAndStick\");
clc; clear all; close all;

%% load data
[Avox,qhat,TE,bvals] = LoadISBI2015Data(1);

%%
[real_params, optim_params, success_rate, min_resnorm] = FitBallAndStick(Avox,qhat,bvals);

model_res = ComputeBallStick(real_params,bvals,qhat)';
SSD = sum((model_res-Avox).^2);

%% compare given values with model values
figure;
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
xlabel('Sample num.')
ylabel('Signal Value')
legend('Data','Model')

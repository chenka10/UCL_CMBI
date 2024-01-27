clear all; clc;

%%
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% perform basic ball and stick fitting

% original voxel
Avox = dwis(:,92,65,72);

% tested voxels
% Avox = dwis(:,28,61,72);
% Avox = dwis(:,82,90,72);
% Avox = dwis(:,55,100,72);

% number of perturbations
N = 10000;

startx = [3.5e+00 3e-03 2.5e-01 pi/2 0];

% setup random noise range to fit parameter values
S0_range = 3e3;
d_range = 0.2;
f_range = 0.5;
theta_range = pi/2;
phi_range = pi;
noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

% perform N ball and stick fitting with random perturbations
[starting_values,fitted_params,resnorms] = RandomBallStickFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store params for best fit
[S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(fitted_params(min_resnorm_index,:));

% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<1)/N;
disp(['success rate: ' num2str(success_rate)]);

%% compute voxel values based on basic model
parameter_hat = fitted_params(min_resnorm_index,:);

model_res = ComputeBallStick_Constrained(parameter_hat,bvals,qhat);

%% compare given values with model values
figure;
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
legend('Data','Model')

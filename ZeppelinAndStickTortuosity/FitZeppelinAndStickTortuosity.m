function [real_params, optim_params, success_rate, min_resnorm] = FitZeppelinAndStickTortuosity(Avox,qhat,bvals)

% number of perturbations
N = 50;

startx = [3.5e+00 3e-03 2.5e-01 0 0];

% setup random noise range to fit parameter values
S0_range = 100;
lam1_range = 100;
f_range = 5;
theta_range = (pi/2)*10;
phi_range = pi*10;
noise_range = [S0_range, lam1_range, f_range, theta_range, phi_range];

% perform N ball and stick fitting with random perturbations
[~,fitted_params,resnorms] = RandomZeppelinStickTortuosityFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
% disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store params for best fit
optim_params = fitted_params(min_resnorm_index,:);
[S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(optim_params);
real_params = [S0,diff,f,theta,phi];

% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<1)/N;
% disp(['success rate: ' num2str(success_rate)]);
end
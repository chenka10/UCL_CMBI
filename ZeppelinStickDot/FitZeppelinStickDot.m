function [real_params, optim_params, success_rate, min_resnorm] = FitZeppelinStickDot(Avox,qhat,bvals)

% number of perturbations
N = 20;

startx = [3.5e+00 3e-03 3e-03 2.5e-01 2.5e-01 pi/2 0];

% setup random noise range to fit parameter values
S0_range = 100;
lam1_range = 100;
lam2_range = 100;
f1_range = 5;
f2_range = 5;
theta_range = (pi/2)*10;
phi_range = pi*10;
noise_range = [S0_range, lam1_range, lam2_range, f1_range, f2_range, theta_range, phi_range];

% perform N ball and stick fitting with random perturbations
[~,fitted_params,resnorms] = RandomZeppelinStickDotFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
% disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store params for best fit
optim_params = fitted_params(min_resnorm_index,:);
real_params = GetRealParamsFromOptimParams_ZeppelinStickDot(fitted_params(min_resnorm_index,:));

% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<1)/N;
% disp(['success rate: ' num2str(success_rate)]);
end
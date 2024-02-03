function [real_params, optim_params, success_rate, min_resnorm] = FitBallAndStick(Avox,qhat,bvals)

% number of perturbations
N = 100;

startx = [3.5e+00 3e-03 2.5e-01 0 0];

% setup random noise range to fit parameter values
S0_range = 10;
d_range = 10;
f_range = 0.5;
theta_range = pi;
phi_range = pi;
noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

% perform N ball and stick fitting with random perturbations
[~,fitted_params,resnorms,~] = RandomBallStickFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store params for best fit
optim_params = fitted_params(min_resnorm_index,:);
[S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(optim_params);
real_params = [S0,diff,f,theta,phi];

% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<0.001)/N;
disp(['success rate: ' num2str(success_rate,5)]);
end
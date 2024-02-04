function [real_params, optim_params, success_rate, min_resnorm] = FitTensorStickDot(Avox,qhat,bvals)

% number of perturbations
N = 150;

startx = [3.5e+00 3e-03 3e-03 3e-03 3e-03 2.5e-01 2.5e-01 0 0 0];

% setup random noise range to fit parameter values
S0_range = 10;
d_s_range = 10;
d_p_range = 10;
d_1_range = 10;
d_2_range = 10;
f1_range = 5;
f2_range = 5;
alpha_range = pi;
beta_range = pi;
gamma_range = pi;
noise_range = [S0_range, d_s_range, d_p_range, d_1_range,d_2_range, f1_range, f2_range, alpha_range, beta_range, gamma_range];

% perform N ball and stick fitting with random perturbations
[~,fitted_params,resnorms] = RandomTensorStickDotFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
% disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store params for best fit
optim_params = fitted_params(min_resnorm_index,:);
real_params = GetRealParamsFromOptimParams_TensorStickDot(fitted_params(min_resnorm_index,:));

% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<1)/N;
% disp(['success rate: ' num2str(success_rate)]);
end
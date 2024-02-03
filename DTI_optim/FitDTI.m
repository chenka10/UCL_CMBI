function [real_params, optim_params, success_rate, min_resnorm] = FitDTI(Avox,qhat,bvals)

% number of perturbations
N = 20;

startx = [3.5e+00 0 0 0 0 0 0];

% setup random noise range to fit parameter values
S0_range = 1;
tensor_range = 0.00001;
noise_range = [S0_range, tensor_range,tensor_range,tensor_range,tensor_range,tensor_range,tensor_range];

% perform N ball and stick fitting with random perturbations
[~,fitted_params,resnorms] = RandomDtiFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
% disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<1)/N;
% disp(['success rate: ' num2str(success_rate)]);

real_params = fitted_params(min_resnorm_index,:);
optim_params = real_params;
end
function [real_params, optim_params, success_rate, min_resnorm] = FitBallAndStickT2(Avox,qhat,bvals,TE,start_params,num_tries)

% number of perturbations
N = 500;

startx = [3.5e+00 3e-03 2.5e-01 0 0];

if nargin>=5
    startx(1:5) = start_params;
end

if nargin>=6
    N = num_tries;
end


% setup random noise range to fit parameter values
S0_range = 100;
d_range = 0.01;
f_range = 0.5;
theta_range = pi;
phi_range = pi;
noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

% perform N ball and stick fitting with random perturbations
[~,fitted_params,resnorms,~] = RandomBallStickT2Fitting(startx,noise_range,Avox,qhat,bvals,TE,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
% disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store params for best fit
optim_params = fitted_params(min_resnorm_index,:);
real_params = GetRealParamsFromOptimParams_BallAndStickT2(optim_params);


% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<0.001)/N;
% disp(['success rate: ' num2str(success_rate,5)]);
end
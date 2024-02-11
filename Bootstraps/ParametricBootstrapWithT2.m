function [results,std_results] = ParametricBootstrapWithT2(model_res,qhat,bvals,TE,N_bootstrap_iterations,N_random_perturbations,sigma_residuals,start_params)
results = zeros(N_bootstrap_iterations,6);

% stds of aggregated results of S0, diff and f
std_results = zeros(N_bootstrap_iterations,3);

startx = [3.5e+00 3e-03 2.5e-01 0 0];

if nargin>=8
    startx = start_params;
end


for i=1:N_bootstrap_iterations

    Avox_sample = model_res + randn(size(model_res))*sigma_residuals;    

    % setup random noise range to fit parameter values
    S0_range = 10;
    d_range = 10;
    f_range = 0.5;
    theta_range = pi;
    phi_range = pi;
    noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

    % perform N ball and stick fitting with random perturbations
    [~,fitted_params,resnorms,~] = RandomBallStickT2Fitting(startx,noise_range,Avox_sample,qhat,bvals,TE,N_random_perturbations);

    % store min resnorm
    [min_resnorm, min_resnorm_index] = min(resnorms);

    success_rate = sum(abs(resnorms-min_resnorm)<0.001)/N_bootstrap_iterations;

    % store params for best fit
    [S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(fitted_params(min_resnorm_index,:));

    results(i,:) = [S0,diff,f,theta,phi,min_resnorm];

    std_results(i,:) = [std(results(1:i,1)), std(results(1:i,2)), std(results(1:i,3))];
end
end


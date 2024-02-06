function [results,std_results] = ParametricBootstrap(model_res,qhat,bvals,N,K,sigma_residuals)
results = zeros(N,6);

% stds of aggregated results of S0, diff and f
std_results = zeros(N,3);


for i=1:N

    Avox_sample = model_res + randn(size(model_res))*sigma_residuals;

    startx = [3.5e+00 3e-03 2.5e-01 0 0];

    % setup random noise range to fit parameter values
    S0_range = 5e3;
    d_range = 10;
    f_range = 0.5;
    theta_range = pi;
    phi_range = pi;
    noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

    % perform N ball and stick fitting with random perturbations
    [~,fitted_params,resnorms,~] = RandomBallStickFitting(startx,noise_range,Avox_sample,qhat,bvals,K);

    % store min resnorm
    [min_resnorm, min_resnorm_index] = min(resnorms);

    % store params for best fit
    [S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(fitted_params(min_resnorm_index,:));

    results(i,:) = [S0,diff,f,theta,phi,min_resnorm];

    std_results(i,:) = [std(results(1:i,1)), std(results(1:i,2)), std(results(1:i,3))];
end
end


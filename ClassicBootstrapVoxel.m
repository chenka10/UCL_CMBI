function [results, std_results, success_rates] = ClassicBootstrapVoxel(Avox,qhat,bvals, Y_pinv,N,K)

x = Y_pinv*log(Avox);

D_DTI = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];
S0_DTI = exp(x(1));

R = D_DTI/trace(D_DTI);
FA_DTI = sqrt(0.5*(3-1/trace(R^2)));

% number of data samples
data_size = size(Avox,1);

results = zeros(N,6);

% stds of aggregated results of S0, diff and f
std_results = zeros(N,3);

success_rates = zeros(N,1);

for i=1:N

    sample_indices = ceil(rand(1,data_size)*data_size);
    Avox_sample = Avox(sample_indices);
    qhat_sample = qhat(:,sample_indices);
    bvals_sample = bvals(sample_indices);

    startx = [S0_DTI trace(D_DTI)/3 FA_DTI 0 0];

   % setup random noise range to fit parameter values
    S0_range = 100;
    d_range = 0.1;
    f_range = 0.1;
    theta_range = pi;
    phi_range = pi;
    noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

    % perform N ball and stick fitting with random perturbations
    [~,fitted_params,resnorms,~] = RandomBallStickFitting(startx,noise_range,Avox_sample,qhat_sample,bvals_sample,K);

    % store min resnorm
    [min_resnorm, min_resnorm_index] = min(resnorms);
    % disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

    % store params for best fit
    [S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(fitted_params(min_resnorm_index,:));

    success_rates(i) = sum(abs(resnorms-min_resnorm)<5)/N;

    results(i,:) = [S0,diff,f,theta,phi,min_resnorm];

    std_results(i,:) = [std(results(1:i,1)), std(results(1:i,2)), std(results(1:i,3))];
end
end

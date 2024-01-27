clear all; clc;

%% load data
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% perform basic ball and stick fitting

% number of bootstrap iterations
N = 200;

% number of random perturbations
K = 10; % in section 1.1 we found that 10 get us 95% chances of finding the minimum

Avox = dwis(:,92,65,72);

% number of data samples
data_size = size(Avox,1);

results = zeros(N,6);

% stds of aggregated results of S0, diff and f
std_results = zeros(N,3);

for i=1:N

    sample_indices = ceil(rand(1,data_size)*data_size);
    Avox_sample = Avox(sample_indices);
    qhat_sample = qhat(:,sample_indices);
    bvals_sample = bvals(sample_indices);

    % Avox_sample = Avox;
    % qhat_sample = qhat;
    % bvals_sample = bvals;

    startx = [3.5e+00 3e-03 2.5e-01 0 0];

    % setup random noise range to fit parameter values
    S0_range = 3e3;
    d_range = 0.2;
    f_range = 0.5;
    theta_range = pi/2;
    phi_range = pi;
    noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

    % perform N ball and stick fitting with random perturbations
    [starting_values,fitted_params,resnorms] = RandomBallStickFitting(startx,noise_range,Avox_sample,qhat_sample,bvals_sample,K);

    % store min resnorm
    [min_resnorm, min_resnorm_index] = min(resnorms);
    % disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

    % store params for best fit
    [S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(fitted_params(min_resnorm_index,:));

    results(i,:) = [S0,diff,f,theta,phi,min_resnorm];

    std_results(i,:) = [std(results(1:i,1)), std(results(1:i,2)), std(results(1:i,3))];
end


%% Display Standard Deviations
figure();
subplot(1,3,1)
plot(1:N,std_results(:,1));
title('S0 standard deviation')
xlabel('num. bootstrap iterations')

subplot(1,3,2)
plot(1:N,std_results(:,2));
title('diff standard deviation')
xlabel('num. bootstrap iterations')

subplot(1,3,3)
plot(1:N,std_results(:,3));
title('f standard deviation')
xlabel('num. bootstrap iterations')

%% Display Histograms
format short;

param_names = {'S0','diff','f'};

figure;
for i=1:3
    
    sigma = std(results(:,i));
    mu = mean(results(:,i));
    per_low = prctile(results(:,i),2.5);
    per_high = prctile(results(:,i),97.5);
    

    subplot(1,3,i)
    histogram(results(:,i),15, 'Normalization', 'probability');
    title({['Histogram of ' param_names{i}], ...
        ['2 sigma: [' num2str(mu-2*sigma) ', ' num2str(mu+2*sigma) '] (Blue)'], ...
        ['95%: [' num2str(per_low) ', ' num2str(per_high) '] (Red)']});
    xlabel(param_names{i})
    hold on;
    line([per_low, per_low], ylim, 'Color', 'r', 'LineWidth', 2);
    line([per_high, per_high], ylim, 'Color', 'r', 'LineWidth', 2);
    line([mu-2*sigma, mu-2*sigma], ylim, 'Color', 'b', 'LineWidth', 2);
    line([mu+2*sigma, mu+2*sigma], ylim, 'Color', 'b', 'LineWidth', 2);
end


%% compute voxel values based on basic model
model_res = ComputeBallStick_Constrained(parameter_hat,bvals,qhat);

%% compare given values with model values
figure;
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
legend('Data','Model')

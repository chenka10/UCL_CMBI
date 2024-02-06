clear all; clc;

%%
addpath("Bootstraps\");

%% load data
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% compute original residuals
selected_slice = 72;
selected_i = 92;
selected_j = 65;

Avox = dwis(:,selected_i,selected_j,selected_slice);

N = 100;

startx = [3.5e+00 3e-03 2.5e-01 0 0];

% setup random noise range to fit parameter values
S0_range = 5e3;
d_range = 10;
f_range = 0.5;
theta_range = pi;
phi_range = pi;
noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

% perform N ball and stick fitting with random perturbations
[starting_values,fitted_params,resnorms,~] = RandomBallStickFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store params for best fit
parameter_hat = fitted_params(min_resnorm_index,:);
[S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(parameter_hat);

model_res = ComputeBallStick_Constrained(parameter_hat,bvals,qhat)';

residuals = Avox - model_res;

% estimate sigma residuals
K = 108;
N = 5;
sigma_residuals = sqrt(sum(residuals.^2)/(K-N));


%% perform basic ball and stick fitting

% number of bootstrap iterations
N = 300;

% number of random perturbations
K = 30; % in section 1.1 we found that 10 get us 95% chances of finding the minimum

% number of data samples
data_size = size(Avox,1);

names = {'Classic','Parametric','Residual','Wild'};

results = cell(4,1);
std_results = cell(4,1);

[curr_results,curr_std_results] = ClassicBootstrap(Avox,qhat,bvals,N,K);
results{1} = curr_results;
std_results{1} = curr_std_results;
fprintf('done classic\n')

[curr_results,curr_std_results] = ParametricBootstrap(model_res,qhat,bvals,N,K,sigma_residuals);
results{2} = curr_results;
std_results{2} = curr_std_results;
fprintf('done parametric\n')

[curr_results,curr_std_results] = ResidualBootstrap(model_res,qhat,bvals,N,K,residuals);
results{3} = curr_results;
std_results{3} = curr_std_results;
fprintf('done residual\n')

[curr_results,curr_std_results] = WildBootstrap(model_res,qhat,bvals,N,K,residuals);
results{4} = curr_results;
std_results{4} = curr_std_results;
fprintf('done wild\n')

%% Display Histograms
format short;

param_names = {'S0','diff','f'};

figure('Position',[10 10 1500 300]);
sgtitle(['Residuals bootstrap results for voxel: [' num2str(selected_i) ', ' num2str(selected_j) ', ' num2str(selected_slice) ']'])
for i=1:3

    sigma = std(results(:,i));
    mu = mean(results(:,i));
    per_low = prctile(results(:,i),2.5);
    per_high = prctile(results(:,i),97.5);


    subplot(1,3,i)
    histogram(results(:,i),15, 'Normalization', 'probability');
    title({['Histogram of ' param_names{i}], ...
        ['\mu=' num2str(mu,3) ', \sigma=' num2str(sigma,3)],...
        ['2\sigma: [' num2str(mu-2*sigma,3) ', ' num2str(mu+2*sigma,3) '], length: ' num2str(4*sigma,3) ' (Blue)'], ...
        ['95%: [' num2str(per_low,3) ', ' num2str(per_high,3) '], length: ' num2str(per_high-per_low,3) ' (Red)']});
    xlabel(param_names{i})
    hold on;
    line([per_low, per_low], [0 1], 'Color', 'r', 'LineWidth', 2);
    line([per_high, per_high], [0 1], 'Color', 'r', 'LineWidth', 2);
    line([mu-2*sigma, mu-2*sigma], [0 1], 'Color', 'b', 'LineWidth', 2);
    line([mu+2*sigma, mu+2*sigma], [0 1], 'Color', 'b', 'LineWidth', 2);
    ylim([0 0.25])
end

%% error scatter
param_names = {'S0','diff','f'};

figure('Position',[10 10 1500 300]);
sgtitle(['Various Bootstrap Errors (voxel: 92,65,72); 2\sigma(blue), 95%(red)'])
for i=1:3
    subplot(1,3,i);

    for iter=1:4
        curr_results = results{iter};
        sigma = std(curr_results(:,i));
        mu = mean(curr_results(:,i));
        per_low = prctile(curr_results(:,i),2.5);
        per_high = prctile(curr_results(:,i),97.5);

        errorbar(iter,mu,mu-per_low,per_high-mu,'Color','r');
        hold on;
        errorbar(iter,mu,2*sigma,2*sigma,'Color','b');
        hold on;
        scatter(iter,mu,'k','filled');
        hold on;
        set(gca,'xtick',[1:4],'xticklabel',names)
        xlim([0 5]);        
    end

    title(param_names{i});
end


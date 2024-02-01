clear all; clc;

%% load data
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% perform basic ball and stick fitting

selected_slice = 72;
selected_i = 92;
selected_j = 65;

% other tested voxels
% selected_i = 28;
% selected_j = 61;
% selected_i = 82;
% selected_j = 90;
% selected_i = 55;
% selected_j = 100;


% number of MCMC iterations (after stabilization!)
N = 100000;

% samples interval
I = 100;

% interval between storing iterations
stabilization_iterations = 30000;

N = N+stabilization_iterations;

% we assume standard deviation of 200
sigma = 200;

Avox = dwis(:,selected_i,selected_j,selected_slice);

% number of data samples
data_size = size(Avox,1);

results = zeros(N,6);

% acceptance_log
acceptangce_log = zeros(1,N);


% setup random noise range to fit parameter values
S0_range = 30;
d_range = 0.00001;
f_range = 0.003;
theta_range = pi/100;
phi_range = pi/100;
noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

% setup starting params
current_params = [3.5e+00 3e-03 2.5e-01 0 0];

current_model_signals = ComputeBallStick(current_params,bvals,qhat)';
current_log_likelihood = ComputeLogLikelihood(sigma,Avox,current_model_signals);

for i=1:N

    % perturbe params
    noise = randn(1,5).*noise_range;
    params = current_params + noise;

    % limit perturbed params
    params(1) = abs(params(1));
    params(2) = abs(params(2));
    params(3) = mod(params(3),1);

    current_model_signals = ComputeBallStick(params,bvals,qhat)';

    resnorm = sum((Avox - current_model_signals).^2);
    log_likelihood = ComputeLogLikelihood(sigma,Avox,current_model_signals);

    if (exp(log_likelihood-current_log_likelihood) > rand())
        current_params = params;
        current_log_likelihood = log_likelihood;
        acceptangce_log(i) = 1;        
    end
    
    results(i,:) = [current_params,resnorm];

    if mod(i,10000)==0
        disp([num2str(i) '/' num2str(N)])
    end

end

results = results(stabilization_iterations:I:end,:);
acceptance_log_all = mean(acceptangce_log(stabilization_iterations:end));


%% Display parameters

param_names = {'S0','diff','f'};

figure('Position',[100 100 1400 200]);
for i=1:numel(param_names)
    subplot(1,numel(param_names),i)
    plot(results(:,i));
    title(param_names{i})
    xlabel('num. iterations')
end

disp(['min err: ' num2str(min(results(:,6)),2)]);


%% Display Histograms
format short;

param_names = {'S0','diff','f'};

figure('Position',[10 10 1500 300]);
sgtitle(['MCMC results for voxel: [' num2str(selected_i) ', ' num2str(selected_j) ', ' num2str(selected_slice) ']'])
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


%% compute voxel values based on basic model
model_res = ComputeBallStick_Constrained(parameter_hat,bvals,qhat);

%% compare given values with model values
figure;
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
legend('Data','Model')

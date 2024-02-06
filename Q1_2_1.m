clear all; clc;

%%

addpath("BallAndStick\");
addpath("Bootstraps\");

%% load data
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% perform basic ball and stick fitting

selected_slice = 72;

i_s = [92 28 82 55];
j_s = [65 61 90 100];

% number of bootstrap iterations
N = 100;

% number of random perturbations
K = 30; % in section 1.1 we found that 10 get us 95% chances of finding the minimum


results = {zeros(N,6),zeros(N,6),zeros(N,6),zeros(N,6)};

% stds of aggregated results of S0, diff and f
std_results = {zeros(N,3),zeros(N,3),zeros(N,3),zeros(N,3)};

for iter=1:4

    selected_i = i_s(iter);
    selected_j = j_s(iter);

    Avox = dwis(:,selected_i,selected_j,selected_slice);

    [curr_results,curr_std_results] = ClassicBootstrap(Avox,qhat,bvals,N,K);

    results{iter} = curr_results;
    std_results{iter} = curr_std_results;

    fprintf('done i=%d, j=%d\n',selected_i,selected_j);
end


%% Display Standard Deviations
iter=4;
figure();
subplot(1,3,1)
plot(1:N,std_results{iter}(:,1));
title('S0 standard deviation')
xlabel('num. bootstrap iterations')

subplot(1,3,2)
plot(1:N,std_results{iter}(:,2));
title('diff standard deviation')
xlabel('num. bootstrap iterations')

subplot(1,3,3)
plot(1:N,std_results{iter}(:,3));
title('f standard deviation')
xlabel('num. bootstrap iterations')

%% Display Histograms
iter = 1;
format short;

selected_i = i_s(iter);
selected_j = j_s(iter);

curr_results = results{iter};
param_names = {'S0','diff','f'};

figure('Position',[10 10 1500 300]);
sgtitle(['Bootstrap results for voxel: [' num2str(selected_i) ', ' num2str(selected_j) ', ' num2str(selected_slice) ']'])
for i=1:3

    sigma = std(curr_results(:,i));
    mu = mean(curr_results(:,i));
    per_low = prctile(curr_results(:,i),2.5);
    per_high = prctile(curr_results(:,i),97.5);


    subplot(1,3,i)
    histogram(curr_results(:,i),15, 'Normalization', 'probability');
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
    ylim([0 0.32])
end

%% error bars
param_names = {'S0','diff','f'};

figure('Position',[10 10 1500 300]);
sgtitle(['Bootstrap Error Bars All Voxels; 2\sigma(blue), 95%(red)'])
for i=1:3
    subplot(1,3,i);
    
    for iter=1:4
        curr_results = results{iter};
        sigma = std(curr_results(:,i));
        mu = mean(curr_results(:,i));
        per_low = prctile(curr_results(:,i),2.5);
        per_high = prctile(curr_results(:,i),97.5);

        bar(iter,mu);
        hold on;
        errorbar(iter,mu,mu-per_low,per_high-mu,'Color','r');
        hold on;        
        errorbar(iter,mu,2*sigma,2*sigma,'Color','b');
        hold on;
        xticks(1:4);
        xlim([0 5]);
        xlabel('Voxel Number');
        if i==1
            ylim([3000 8000]);
        end
    end
    
    title(param_names{i});
end


%% compute voxel values based on basic model
model_res = ComputeBallStick_Constrained(parameter_hat,bvals,qhat);

%% compare given values with model values
figure;
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
legend('Data','Model')

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


% number of bootstrap iterations
N = 100;

% number of random perturbations
K = 15; % in section 1.1 we found that 10 get us 95% chances of finding the minimum

Avox = dwis(:,selected_i,selected_j,selected_slice);

Y = GetDesignMatrix(qhat,bvals);
Y_pinv = pinv(Y);

tic
[results,std_results,success_rates] = ClassicBootstrapVoxel(Avox,qhat,bvals,Y_pinv,N,K);
toc


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

figure('Position',[10 10 1500 300]);
sgtitle(['Bootstrap results for voxel: [' num2str(selected_i) ', ' num2str(selected_j) ', ' num2str(selected_slice) ']'])
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



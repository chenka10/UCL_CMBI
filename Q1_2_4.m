clear all; clc;

%% load data
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% perform basic ball and stick fitting
downsample_rate = 1;
min_dwis = dwis(:,1:downsample_rate:145,1:downsample_rate:174,:);

size(min_dwis)
%%
figure;
imshow(flipud(squeeze(min_dwis(1,:,:,72)>2000)'), []);


%%
fprintf('voxel num = %d\n',sum(min_dwis(1,:,:,72)>2000,'all'))

N=10000;
stabilization = 500;
I = 10;

size_x = size(min_dwis,2);
size_y = size(min_dwis,3);

results = zeros(size_x,size_y,N+stabilization,6);
acceptance_counts = zeros(size_x,size_y);

% acceptance_log
acceptangce_count = zeros(size_x,size_y);

Y = GetDesignMatrix(qhat,bvals);
Y_pinv = pinv(Y);

selected_slice = 72;
for selected_i=1:size_x
    for selected_j=1:size_y        

        Avox = min_dwis(:,selected_i,selected_j,selected_slice);

        if Avox(1)<2000
            continue
        end
        
        Y = GetDesignMatrix(qhat,bvals);
        Y_pinv = pinv(Y);

        % MCMC(Avox,qhat,bvals,Y_pinv,N)
        tic
        [curr_results,acceptance_count]=MCMC(Avox,qhat,bvals,Y_pinv,N,stabilization);
        toc

        results(selected_i,selected_j,:,:) = curr_results;
        acceptance_counts(selected_i,selected_j) = acceptance_count;

        fprintf('i=%d, j=%d\n',selected_i,selected_j);
    end
end

%% Display parameters

param_names = {'S0','diff','f'};

figure('Position',[100 100 1400 200]);
for i=1:numel(param_names)
    subplot(1,numel(param_names),i)
    plot(squeeze(results(92,65,:,i)));
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

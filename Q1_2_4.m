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

        fprintf('i=%d, j=%d, acc =%d\n',selected_i,selected_j,acceptance_count/N);
    end
end

%% Display Uncertainty results

param_names = {'S0','d','f'};
ranges = [1 0.6 0.6];

figure('Position',[100 100 800 800]);
sgtitle('MCMC Uncertainty Results')
for i=1:numel(param_names)
    subplot(2,2,i)
    img = flipud(squeeze(4*std(results(:,:,stabilization:I:(N+stabilization),i),0,3))');
    imshow(img,[0 ranges(i)*max(img,[],'all')]);
    title(['2\sigma range: ' param_names{i}])  
    colorbar()
end
subplot(2,2,4)
imshow(flipud(acceptance_counts'/N));
colorbar()
title('acceptance rate');

%% process directional historgrams for display

thetas = results(:,:,stabilization:I:(N+stabilization),4);
phis = results(:,:,stabilization:I:(N+stabilization),5);

thetas_mean = squeeze(mean(thetas,3));
phis_mean = squeeze(mean(phis,3));

thetas_res = thetas - thetas_mean;
phis_res = phis - phis_mean;

range_x = (1:145);
range_y = (1:174);

nbins = 7;

phi_bins = zeros(size_x,size_y,nbins);
phi_counts = zeros(size_x,size_y,nbins);

for i=1:size_x
    for j=1:size_y        
        curr_phis_res = squeeze(phis_res(i,j,:));

        if all(curr_phis_res==0)
            continue
        end
        
        [phis_res_binned,curr_phi_bins] = histcounts(curr_phis_res,nbins);

        phi_bins(i,j,:) = (curr_phi_bins(1:nbins)+curr_phi_bins(2:(nbins+1)))/2;
        phi_counts(i,j,:) = phis_res_binned;       

        fprintf('%d,%d\n',i,j);
    end
end

%% plot directional uncertainty

centers = [size_x/2 size_y/2; 73 80; 33 96];
window_sizes = [size_x/2 size_y/2; 10 12; 10 12];

figure('Position',[100 100 1300 500]);
sgtitle('Directional Uncertainty')
phi_counts_reshaped = phi_counts.^0.3;
phi_counts_reshaped = phi_counts_reshaped./sum(phi_counts_reshaped,3);

for fc=1:3
    subplot(1,3,fc);

    range_x = int32((centers(fc,1)-window_sizes(fc,1)+1):(centers(fc,1)+window_sizes(fc,1)));
    range_y = int32((centers(fc,2)-window_sizes(fc,2)+1):(centers(fc,2)+window_sizes(fc,2)));

    for bin=1:nbins
        phi_counts_bin = squeeze(phi_counts_reshaped(range_x,range_y,bin));
        phi_bins_bin = squeeze(phi_bins(range_x,range_y,bin));        

        phis = phis_mean(range_x,range_y) + phi_bins_bin;

        dir_vecs_x = (cos(phis)).*phi_counts_bin*0.5;
        dir_vecs_y = (sin(phis)).*phi_counts_bin*0.5;

        quiver(range_x,range_y,(dir_vecs_x)',(dir_vecs_y)','ShowArrowHead','off','Color','b');
        hold on
        quiver(range_x,range_y,(-dir_vecs_x)',(-dir_vecs_y)','ShowArrowHead','off','Color','b');
        xlim([range_x(1) range_x(end)])
        ylim([range_y(1) range_y(end)])
    end
    if fc==1
        hold on;
        rectangle('EdgeColor','g','LineWidth',2,'Position',[centers(2,1)+1-window_sizes(2,1) centers(2,2)+1-window_sizes(2,2) 2*window_sizes(2,1) 2*window_sizes(2,2)])
        rectangle('EdgeColor','r','LineWidth',2,'Position',[centers(3,1)+1-window_sizes(3,1) centers(3,2)+1-window_sizes(3,2) 2*window_sizes(3,1) 2*window_sizes(3,2)])
    end
    if fc ==2
        hold on;
        rectangle('EdgeColor','g','LineWidth',2,'Position',[centers(2,1)+1-window_sizes(2,1) centers(2,2)+1-window_sizes(2,2) 2*window_sizes(2,1)-1 2*window_sizes(2,2)-1])
    end
    if fc ==3
        hold on;
        rectangle('EdgeColor','r','LineWidth',2,'Position',[centers(3,1)+1-window_sizes(3,1) centers(3,2)+1-window_sizes(3,2) 2*window_sizes(3,1)-1 2*window_sizes(3,2)-1])
    end
    daspect([1 1 1])
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

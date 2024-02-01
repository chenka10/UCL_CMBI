clear all; clc;

%% load data
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% compute maps

res_map = zeros(size(dwis,1),size(dwis,2),6);
N = 30;
selected_slice = 72;

for i=1:145
    for j=1:174
        try
            % load values for chosen voxel
            Avox = dwis(:,i,j,72);

            % perform basic ball and stick fitting
            startx = [3.5e+00 3e-03 2.5e-01 0 0];

            % setup random noise range to fit parameter values
            S0_range = 5e3;
            d_range = 10;
            f_range = 0.5;
            theta_range = pi;
            phi_range = pi;
            noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

            % fit ball and stick using random perturbations
            [starting_values,fitted_params,resnorms,~] = RandomBallStickFitting(startx,noise_range,Avox,qhat,bvals,N);

            % find minimal resnorm
            [min_resnorm, min_resnorm_index] = min(resnorms);            

            % compute voxel values based on basic model
            parameter_hat = fitted_params(min_resnorm_index,:);

            % Extract the parameters
            [S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(parameter_hat);

            res_map(i,j,:) = [S0,diff,f,theta,phi,min_resnorm];            
            
            if (mod(j,50)==1)
                disp(['i: ' num2str(i) ',j: ' num2str(j)])
            end
        catch exception            
            fprintf('Error message: %s\n', exception.message);
        end
    end
end

%%
save('results_30_new.mat','res_map');

%% Displaying maps
figure('Position',[100 100 550 500]);
sgtitle(['Slice: ' num2str(selected_slice)])
subplot(2,2,1)
imshow(flipud(res_map(:,:,1)'), []);
xticks([])
yticks([])
colorbar()
title('S0 Map')
subplot(2,2,2)
imshow(flipud(res_map(:,:,2)'), [0 0.005]);
xticks([])
yticks([])
colorbar()
title('diff Map')
subplot(2,2,3)
imshow(flipud(res_map(:,:,3)'), []);
xticks([])
yticks([])
colorbar()
title('f Map')
subplot(2,2,4)
imshow(flipud(res_map(:,:,6)'), [0 2e7]);
xticks([])
yticks([])
colorbar()
title('Residual Error Map')

%% Displaying diraction map

thetas = res_map(:,:,4);
phis= res_map(:,:,5);

dir_vecs_x = (sin(thetas).*cos(phis)).*res_map(:,:,3);
dir_vecs_y = (sin(thetas).*sin(phis)).*res_map(:,:,3);

figure;
subplot(1,3,1)
quiver(1:145,1:174,(dir_vecs_x)',flipud(dir_vecs_y)','ShowArrowHead','off');
daspect([1 1 1])

subplot(1,3,2)
quiver(1:145,1:174,(dir_vecs_x)',flipud(dir_vecs_y)','ShowArrowHead','off');
daspect([1 1 1])

subplot(1,3,3)
quiver(1:145,1:174,(dir_vecs_x)',flipud(dir_vecs_y)','ShowArrowHead','off');
daspect([1 1 1])



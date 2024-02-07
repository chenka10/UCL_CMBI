clear all; clc;

%%

addpath("BallAndStick\");
addpath("DTI_optim\");
addpath("ZeppelinAndStick\");
addpath("ZeppelinAndStickTortuosity\");
addpath("ZeppelinStickDot\");
addpath("ZeppelinTwoSticks\");
addpath("TensorStickDot\");

%% load data
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);



% downsample dwis
downsample_step = 1;
dwis = dwis(:,1:downsample_step:size(dwis,2),1:downsample_step:size(dwis,3),:);

%% compute maps

% params for AIC
Ns = [6 7 6];

% models fitting functions
fit_funcs = {@FitBallAndStick,@FitZeppelinAndStick,@FitZeppelinAndStickTortuosity};
M = numel(fit_funcs);

res_params = {zeros(size(dwis,2),size(dwis,3),7),...
    zeros(size(dwis,2),size(dwis,3),8),...
    zeros(size(dwis,2),size(dwis,3),7)};

res_map = zeros(size(dwis,2),size(dwis,3));
res_map_values = zeros(size(dwis,2),size(dwis,3),M);
selected_slice = 72;

Y = GetDesignMatrix(qhat,bvals);
Y_pinv = pinv(Y);

num_tries = [20 15 15];

for j=1:size(dwis,3)
    for i=1:size(dwis,2)    
        tic
        Avox = dwis(:,i,j,selected_slice);
        if max(Avox)<2000
            continue
        end

        x = Y_pinv*log(Avox);

        D_DTI = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];
        S0_DTI = exp(x(1));

        R = D_DTI/trace(D_DTI);
        FA_DTI = sqrt(0.5*(3-1/trace(R^2)));

        dti_results = [S0_DTI,trace(D_DTI)/3, FA_DTI];

        best_AIC = 99999;        

        for m=1:M
            fit_func = fit_funcs{m};
            [real_params, optim_params, success_rate, min_resnorm] = fit_func(Avox,qhat,bvals,dti_results,num_tries(m));

            K = numel(Avox);
            N = Ns(m);

            fprintf('m=%d, success=%s\n',m,num2str(success_rate,4));

            AIC = ComputeAIC(N,K,min_resnorm) + (2*N*(N+1))/(K-N-1);

            res_map_values(i,j,m) = AIC;

            if AIC<best_AIC
                best_AIC = AIC;
                res_map(i,j) = m;
            end

            res_params{m}(i,j,:) = [real_params success_rate min_resnorm];
        end

        toc
        fprintf('i=%d, j=%d, best=%d\n',i,j,res_map(i,j));
    end
end

%%
save('results_30_new.mat','res_map');

%% Displaying maps
selected_slice=72;
res_params_map = res_params{1};
figure('Position',[100 100 550 500]);
sgtitle(['Slice: ' num2str(selected_slice)])
subplot(2,2,1)
imshow(flipud(res_params_map(:,:,1)'), []);
xticks([])
yticks([])
colorbar()
title('S0 Map')
subplot(2,2,2)
imshow(flipud(res_params_map(:,:,2)'), [0 0.005]);
xticks([])
yticks([])
colorbar()
title('diff Map')
subplot(2,2,3)
imshow(flipud(res_params_map(:,:,3)'), []);
xticks([])
yticks([])
colorbar()
title('f Map')
subplot(2,2,4)
imshow(flipud(res_params_map(:,:,end)'), []);
xticks([])
yticks([])
colorbar()
title('Residual Error Map')

%% display model selection

err_img = (flipud(res_map_values(:,:,1)') - flipud(res_map_values(:,:,3)'));
err_img = err_img - err_img(1,1,1);
norm_err_img = NormalizeImg(err_img);
norm_err_img = norm_err_img - norm_err_img(1,1,1);
norm_err_img(norm_err_img>0.1) = 0.1;
norm_err_img(norm_err_img<-0.1) = -0.1;
norm_err_img = NormalizeImg(norm_err_img);

selected_model_map = (flipud(res_map'));

f_map = flipud(res_params{3}(:,:,3)');

figure('Position',[100 100 800 700]);
ax1=subplot(2,2,1);
imshow(selected_model_map,[]);
colormap(ax1,"jet");
xticks([])
yticks([])
title({'Selected Model'})
subtitle('1: BallStick, 2: ZeppelinStick, 3: ZeppelinStock (Tortuosity)')
colorbar()

subplot(2,2,2)
imshow(f_map, []);
xticks([])
yticks([])
title('f coefficient')
colorbar()


h_sp=subplot(2,2,3);
err_img_mean = nanmean(err_img,'all');
range = max(err_img,[],'all') - min(err_img,[],'all');
h = imshow(err_img,[err_img_mean-range*0.1 err_img_mean+range*0.1]);
xticks([])
yticks([])
colorbar()

mod_spring_range=(0:(1/255):1) - 0.5;
mod_spring = spring;
mod_spring = mod_spring.*repmat(abs(mod_spring_range),3,1)'*2;
% mod_spring(mod_spring_range>0.1,:) = 0;

colormap(h_sp,mod_spring)
title({'Difference: AIC(BallStick) - AIC(ZeppelinSitkc(Tortuosity))'})
subtitle('Lower -> BallStick was better')

subplot(2,2,4)
h_f = imshow(f_map,[]);
hold on;
h_err = image(grayscaleToColor(norm_err_img,spring));
n = sqrt(sum((err_img-err_img(1,1,:)).^2,3));
h_err.AlphaData = abs(norm_err_img - norm_err_img(1,1,1))*5;
xticks([])
yticks([])
title('Difference on f')



%%
best_res_index_map = zeros(size(res_map));

for i=1:size(res_map,1)
    for j=1:size(res_map,2)
        curr_min = inf;
        for m=1:3
            curr_model_resnorm = res_params{m}(i,j,end);

            if curr_model_resnorm<curr_min
                curr_min = curr_model_resnorm;
                best_res_index_map(i,j)=m;
            end
        end        
    end
end

%% display model selection

figure;
imshow((flipud(best_res_index_map')),[]);

function img = NormalizeImg(input)
input = input - min(input,[],'all');
img = input./max(input,[],'all');
end

function colorImage = grayscaleToColor(grayImage, colormapName)
    % Convert grayscale image to indexed image
    indexedImage = gray2ind(grayImage, 256);
    
    % Apply colormap to indexed image
    rgbImage = ind2rgb(indexedImage, colormapName);
    
    % Convert to uint8 data type
    colorImage = uint8(rgbImage * 255);
end

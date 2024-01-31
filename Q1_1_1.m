%% loading data
load('data');

dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% generating maps

% params
selected_slice = 72;
img_width = size(dwis,2);
img_height = size(dwis,3);

% making empty maps
meanDiff_map = zeros(img_width,img_height);
FA_map = zeros(img_width,img_height);
dir_map = zeros(img_width,img_height,3);

% compute design matrix
Y = GetDesignMatrix(qhat,bvals);
Y_pinv = pinv(Y);

% iterating on all voxels to fill maps
for i=1:img_width
    for j=1:img_height
        if min(dwis(:,i,j,selected_slice))<=0
            continue
        end
        D = GetDiffusionMatrixForVoxel(dwis(:,i,j,selected_slice),Y_pinv);
        meanDiff_map(i,j) = (trace(D))/3;
        [V,M] = eig(D);

        R = D/trace(D);
        FA_map(i,j) = sqrt(0.5*(3-1/trace(R^2)));

        % Find the index of the maximum eigenvalue
        [~, maxEigenvalueIndex] = max(diag(M));
        
        % Extract the corresponding eigenvector
        mainEigenvector = V(:, maxEigenvalueIndex);

        dir_map(i,j,:) = abs(mainEigenvector)*FA_map(i,j);
    end
end

%% Display S0 slice
figure;
imshow(flipud(squeeze(dwis(1,:,:,72))'), []);

%% Display Maps
figure('Position',[100 100 1000 300]);
sgtitle(['Slice: ' num2str(selected_slice)])
subplot(1,3,1)
imshow(flipud(meanDiff_map'), []);
xticks([])
yticks([])
colorbar()
title('Mean Diffusion')
subplot(1,3,2)
imshow(flipud(FA_map'), []);
xticks([])
yticks([])
colorbar()
title('Fractional Anisotropy')
subplot(1,3,3)
imshow(flipud(permute(dir_map,[2,1,3])), []);
xticks([])
yticks([])
colorbar()
title('Directional Map')


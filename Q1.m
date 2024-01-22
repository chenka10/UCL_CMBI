load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%%

% Middle slice of the 1st image volume, which has b=0
imshow(flipud(squeeze(dwis(1,:,:,72))'), []);

figure;
% Middle slice of the 2nd image volume, which has b=1000
imshow(flipud(squeeze(dwis(2,:,:,72))'), []);


%%
selected_slice = 72;
img_width = size(dwis,2);
img_height = size(dwis,3);

meanDiff_map = zeros(img_width,img_height);
FA_map = zeros(img_width,img_height);
dir_map = zeros(img_width,img_height,3);

Y = GetDesignMatrix(qhat,bvals);
Y_pinv = pinv(Y);

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

%%
figure;
imshow(flipud(meanDiff_map'), []);

%%
figure;
imshow(flipud(FA_map'), []);

%%
figure;
imshow(flipud(permute(dir_map,[2,1,3])), []);
%%
figure;
quiver3(flipud(permute(dir_map,[2,1,3])));

function Y = GetDesignMatrix(qhat,b)
    images_num = size(qhat,2);
    qx = qhat(1,:);
    qy = qhat(2,:);
    qz = qhat(3,:);
    
    Y = [ones(1,images_num); -b.*qx.^2; -2*b.*qx.*qy; -2*b.*qx.*qz; -b.*qy.^2; -2*b.*qy.*qz; -b.*qz.^2]';
end

function D = GetDiffusionMatrixForVoxel(Avox,Y_pinv)

x = Y_pinv*log(Avox);

D = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];
end


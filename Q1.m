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

% Avox = dwis(:,92,65,72);
Avox = dwis(:,73,110,72);

%%
images_num = size(Avox,1);
qx = qhat(1,:);
qy = qhat(2,:);
qz = qhat(3,:);
b = bvals;

Y = [ones(1,images_num); -b.*qx.^2; -2*b.*qx.*qy; -2*b.*qx.*qz; -b.*qy.^2; -2*b.*qy.*qz; -b.*qz.^2]';
%%
x = Y\log(Avox);
%%
D = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];


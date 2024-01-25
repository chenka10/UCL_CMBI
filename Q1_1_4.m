clear all; clc;

%%
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

% load values for chosen voxel
Avox = dwis(:,92,65,72);

%% perform basic ball and stick fitting

N = 1000;

starting_values = zeros(N,5);
fitted_params = zeros(N,5);
resnorms = zeros(N,1);

for i=1:N

startx = [3.5e+00 3e-03 2.5e-01 pi/2 0];

S0_range = 3e3;
d_range = 0.1;
f_range = 0.5;
theta_range = pi/2;
phi_range = pi;

noise = randn(1,5).*[S0_range, d_range, f_range, theta_range, phi_range];
startx = startx + noise;
startx(1) = abs(startx(1));
startx(2) = abs(startx(2));
startx(3) = max(min(startx(3),1),0);
startx(4) = max(min(startx(4),pi),0);

starting_values(i,:) = startx;

startx(1) = sqrt(startx(1));
startx(2) = sqrt(startx(2));
startx(3) = log(startx(3)/(1-startx(3)));
startx(4) = tan(startx(4)-pi/2);

h=optimset('MaxFunEvals',20000,...
 'Algorithm','quasi-newton',...
 'MaxIter',200,...
 'TolX',1e-10,...
 'TolFun',1e-10, ...
 'Display','none');

[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_Constrained',startx,h,Avox,bvals,qhat);

fitted_params(i,:) = parameter_hat;
resnorms(i) = RESNORM;

end

[min_resnorm, min_resnorm_index] = min(resnorms);
disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

%% compute voxel values based on basic model

parameter_hat = fitted_params(min_resnorm_index,:);

model_res = ComputeBallStick_Constrained(parameter_hat,bvals,qhat);

%% compare given values with model values
figure;
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
legend('Data','Model')

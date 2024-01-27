clear all; clc;

%% load data
load('data');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

%% perform basic ball and stick fitting
Avox = dwis(:,92,65,72);

startx = [3.5e+00 3e-03 2.5e-01 pi/2 0];

% inverse stating params
startx = GetOptimParamsFromRealParams(startx);

h=optimset('MaxFunEvals',20000,...
 'Algorithm','quasi-newton',...
 'MaxIter',200,...
 'TolX',1e-10,...
 'TolFun',1e-10, ...
 'Display','final');

[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_Constrained',startx,h,Avox,bvals,qhat);

%% compute voxel values based on basic model
model_res = ComputeBallStick_Constrained(parameter_hat,bvals,qhat);

%% compare given values with model values
figure;
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
legend('Data','Model')

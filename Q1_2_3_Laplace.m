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

% number of random perturbations
K = 30; % in section 1.1 we found that 10 get us 95% chances of finding the minimum

Avox = dwis(:,selected_i,selected_j,selected_slice);

% number of data samples
data_size = size(Avox,1);

startx = [4257,0.001,0.35,0.981,-2.5];

h=optimset('MaxFunEvals',20000,...
    'Algorithm','quasi-newton',...
    'MaxIter',2000,...
    'TolX',1e-10,...
    'TolFun',1e-10, ...
    'Display','none');

[parameter_hat,RESNORM,~,~,~,hessian]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);

disp(['min SSD: ' num2str(RESNORM) ]);

%% compute standard deviations using Hessian
sigma = 200;
real_hessian = (1/(2*sigma^2))*hessian;

std_values = sqrt(diag(inv(real_hessian)));

disp('sigma Values: [S0,diff,f,theta,phi]:')
disp(std_values)
 



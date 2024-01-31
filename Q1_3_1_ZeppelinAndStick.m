%% clearing

clc; clear all; close all;

%% load data

% Load the diffusion signal
fid = fopen('isbi2015_data_normalised.txt', 'r', 'b');
fgetl(fid); % Read in the header
D = fscanf(fid, '%f', [6, inf])'; % Read in the data
fclose(fid);
% Select the first of the 6 voxels
Avox = D(:,1);
% Load the protocol
fid = fopen('isbi2015_protocol.txt', 'r', 'b');
fgetl(fid);
A = fscanf(fid, '%f', [7, inf]);
fclose(fid);
% Create the protocol
qhat = A(1:3,:);
G = A(4,:);
delta = A(5,:);
smalldel = A(6,:);
TE = A(7,:);
GAMMA = 2.675987E8;
bvals = ((GAMMA*smalldel.*G).^2).*(delta-smalldel/3);
% convert bvals units from s/m^2 to s/mm^2
bvals = bvals/10^6;


%% perform basic ball and stick fitting

% number of perturbations
N = 100;

startx = [3.5e+00 3e-03 3e-03 2.5e-01 pi/2 0];

% setup random noise range to fit parameter values
S0_range = 10;
lam1_range = 1;
lam2_range = 1;
f_range = 0.5;
theta_range = pi/2;
phi_range = pi;
noise_range = [S0_range, lam1_range, lam2_range, f_range, theta_range, phi_range];

% perform N ball and stick fitting with random perturbations
[starting_values,fitted_params,resnorms] = RandomZeppelinStickFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store params for best fit
params_real = GetRealParamsFromOptimParams_ZeppelinStick(fitted_params(min_resnorm_index,:));

% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<1)/N;
disp(['success rate: ' num2str(success_rate)]);

%% compute voxel values based on basic model
parameter_hat = fitted_params(min_resnorm_index,:);

model_res = ComputeBallStick_Constrained(parameter_hat,bvals,qhat);

%% compare given values with model values
figure;
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
legend('Data','Model')

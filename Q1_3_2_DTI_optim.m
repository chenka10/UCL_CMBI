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

startx = [3.5e+00 0 0 0 0 0 0];

% setup random noise range to fit parameter values
S0_range = 1;
tensor_range = 0.00001;
noise_range = [S0_range, tensor_range,tensor_range,tensor_range,tensor_range,tensor_range,tensor_range];

% perform N ball and stick fitting with random perturbations
[starting_values,fitted_params,resnorms] = RandomDtiFitting(startx,noise_range,Avox,qhat,bvals,N);

% store min resnorm
[min_resnorm, min_resnorm_index] = min(resnorms);
disp(['min SSD: ' num2str(min_resnorm) ', at iter: ' num2str(min_resnorm_index)]);

% store success rate
success_rate = sum(abs(resnorms-min_resnorm)<1)/N;
disp(['success rate: ' num2str(success_rate)]);

parameter_hat = fitted_params(min_resnorm_index,:);

%% compute voxel values based on basic model
model_res_DTI_optim = ComputeDti(parameter_hat,bvals,qhat);

%% compare given values with model values
figure('Position',[100 100 1500 400]);
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res_DTI_optim, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
xlabel('Sample num.')
ylabel('Signal Value')
title(['DTI Match Plot; SSD=' num2str(min_resnorm)])
legend('Data','Model')

%% clearing
addpath("./ZeppelinAndStick");
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

%%
[real_params, optim_params, success_rate, min_resnorm] = FitZeppelinAndStick(Avox,qhat,bvals);

model_res = ComputeZeppelinStick(real_params,bvals,qhat)';
SSD = sum((model_res-Avox).^2);

%% compare given values with model values
figure('Position',[100 100 1500 400]);
plot(Avox, ' bs', 'MarkerSize', 6, 'LineWidth', 2);
hold on;
plot(model_res, ' rx', 'MarkerSize', 6, 'LineWidth', 2);
xlabel('Sample num.')
ylabel('Signal Value')
title(['Zeppelin and Stick Match Plot; SSD=' num2str(min_resnorm)])
legend('Data','Model')

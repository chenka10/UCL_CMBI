%% clearing
clc; clear all; close all;

%%
addpath("BallAndStickT2\")
addpath("BallAndStick\")
addpath("experimentDesign\")

%% load ISBI2015 data
selected_voxel = 1;
[Avox,qhat,TE,bvals] = LoadISBI2015Data(selected_voxel);

%% compute shells
shells_idencies = zeros(36,121);

unique_bvals = unique(bvals);
unique_bvals = unique_bvals(2:end); % remove b=0
TE_per_shell = zeros(numel(unique_bvals),1);
for i=1:numel(unique_bvals)    
    idencies_range = 1:numel(bvals);
    idencies = idencies_range(bvals==unique_bvals(i));
    shell_TE = TE(idencies(1));
    TE_per_shell(i) = shell_TE;
    a = bvals==unique_bvals(i);
    b = (bvals==0).*(TE==shell_TE);
    all_idencies = idencies_range(logical(a+b));
    shells_idencies(i,:) = all_idencies;   

end

figure
subplot(2,1,1)
plot(unique_bvals);
title('b values')
xlabel('shell number');
ylabel('b [s/mm^2]')
subplot(2,1,2)
plot(TE_per_shell);
title('TE values')
xlabel('shell number');
ylabel('TE [s]')

%% generate shell couples

N = 36; % Define the maximum number in the array
num_elements = N * (N - 1) / 2; % Calculate the number of unique couples

couples = zeros(num_elements, 2); % Initialize array to store couples

idx = 1; % Initialize index for storing couples

for i = 1:N-1
    for j = i+1:N
        couples(idx, :) = [i, j]; % Store the current couple
        idx = idx + 1;
    end
end

shells_idencies_couples = zeros(num_elements,2*size(shells_idencies,2));

for i=1:num_elements
shells_idencies_couples(i,:) = [shells_idencies(couples(i,1),:) shells_idencies(couples(i,2),:)];
end

shells_idencies = shells_idencies_couples;

%% Compute Fisher matrix for Ball and Stick on all data

params = [1.009865516767256   0.001432055891641   0.574928310559141  1.59686173476731723846   6.200322564900709];
params_norm = params(1:3)'*params(1:3);

F = ComputeFisherMatrix(params,bvals,qhat).*params_norm;

%% Compute Fisher matrix for Ball and Stick on all data (with T2)
T2 = 0.07; % [s]
params = [1.009865516767256   0.001432055891641   0.574928310559141  1.59686173476731723846   6.200322564900709 T2];
params_norm = [params(1:3) params(6)]'*[params(1:3) params(6)];

F = ComputeFisherMatrixWithT2(params,bvals,qhat,TE).*params_norm;



%% Find optimal shell 
params = [1.009865516767256   0.001432055891641   0.574928310559141  -1.544730918822476   6.200322564900709];
params_norm = params(1:3)'*params(1:3);

[A_optim,D_optim,E_optim,T_optim,eigenvalues] = FindOptimShell(shells_idencies, params,params_norm,@ComputeFisherMatrix,bvals,qhat,TE);

[A_optim_val,A_optim_val_i] = min(A_optim);
[D_optim_val,D_optim_val_i] = min(D_optim);
[E_optim_val,E_optim_val_i] = max(E_optim);
[T_optim_val,T_optim_val_i] = max(T_optim);

%% Find optimal shell (With T2)
params = [1.009865516767256   0.001432055891641   0.574928310559141  -1.544730918822476   6.200322564900709 0.07];
params_norm = [params(1:3) params(6)]'*[params(1:3) params(6)];

[A_optim,D_optim,E_optim,T_optim,eigenvalues] = FindOptimShell(shells_idencies, params,params_norm,@ComputeFisherMatrixWithT2,bvals,qhat,TE);

[A_optim_val,A_optim_val_i] = min(A_optim);
[D_optim_val,D_optim_val_i] = min(D_optim);
[E_optim_val,E_optim_val_i] = max(E_optim);
[T_optim_val,T_optim_val_i] = max(T_optim);

%% plot optimality values

figure('Position',[100 100 1000 700]);
subplot(2,2,1)
plot(A_optim);
xlabel('shell number')
title({'trace(F^{-1}) - minimize for A-optimality',['min: ' num2str(A_optim_val) ',min_i: ' num2str(A_optim_val_i)]});
subplot(2,2,2)
plot(D_optim);
xlabel('shell number')
title({'det(F^{-1}) - minimize for D-optimality',['min: ' num2str(D_optim_val) ',min_i: ' num2str(D_optim_val_i)]});
subplot(2,2,3)
plot(E_optim);
xlabel('shell number')
title({'min eigenvalue of F - maximize for E-optimality',['max: ' num2str(E_optim_val) ',max_i: ' num2str(E_optim_val_i)]});
subplot(2,2,4)
plot(T_optim);
xlabel('shell number')
title({'trace(F) - maximiaze for T-optimality',['max: ' num2str(T_optim_val) ',max_i: ' num2str(T_optim_val_i)]});

%% plot eigenvalues progression

figure;
for i=1:size(eigenvalues,2)
subplot(2,2,i)
plot(eigenvalues(:,i));
[~,optim_val] = max(eigenvalues(:,i));
xlabel('shell number')
title({['eigenvalue num ' num2str(i)],['max: ' num2str(optim_val)]})
end












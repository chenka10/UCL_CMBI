%% clearing
clc; clear all; close all;

%% load ISBI2015 data
selected_voxel = 1;
[Avox,qhat,TE,bvals] = LoadISBI2015Data(selected_voxel);

%% compute shells
shells_idencies = zeros(36,121);

unique_bvals = unique(bvals);
unique_bvals = unique_bvals(2:end); % remove b=0
for i=1:numel(unique_bvals)    
    idencies_range = 1:numel(bvals);
    idencies = idencies_range(bvals==unique_bvals(i));
    shell_TE = TE(idencies(1));
    a = bvals==unique_bvals(i);
    b = (bvals==0).*(TE==shell_TE);
    all_idencies = idencies_range(logical(a+b));
    shells_idencies(i,:) = all_idencies;   

end

%% Compute Fisher matrix for Ball and Stick on all data

params = [1.009865516767256   0.001432055891641   0.574928310559141  1.59686173476731723846   6.200322564900709];
params_norm = params(1:3)'*params(1:3);

F = ComputeFisherMatrix(params,bvals,qhat).*params_norm;



%% Find optimal shell
params = [1.009865516767256   0.001432055891641   0.574928310559141  -1.544730918822476   6.200322564900709];
params_norm = params(1:3)'*params(1:3);

A_optim = zeros(size(shells_idencies,1),1);
D_optim = zeros(size(shells_idencies,1),1);
E_optim = zeros(size(shells_idencies,1),1);
T_optim = zeros(size(shells_idencies,1),1);

eigenvalues = zeros(size(shells_idencies,1),3);


for i=1:size(shells_idencies,1)
    curr_shell_idencies = shells_idencies(i,:);
    curr_bvals = bvals(curr_shell_idencies);
    curr_qhat = qhat(curr_shell_idencies);
    F = ComputeFisherMatrix(params,curr_bvals,curr_qhat).*params_norm;
    F_inv = pinv(F);
    A_optim(i)= trace(F_inv);
    D_optim(i) = det(F_inv);

    curr_eigenvalues = eig(F);
    min_eigenvalue = min(curr_eigenvalues);
    E_optim(i) = min_eigenvalue;
    
    eigenvalues(i,:) = curr_eigenvalues;

    T_optim(i) = F(1,1);
end

[~,A_optim_val] = min(A_optim);
[~,D_optim_val] = min(D_optim);
[~,E_optim_val] = max(E_optim);
[~,T_optim_val] = max(T_optim);

%% plot optimality values

figure('Position',[100 100 1000 700]);
subplot(2,2,1)
plot(A_optim);
xlabel('shell number')
title({'trace(F^{-1}) - minimize for A-optimality',['min: ' num2str(A_optim_val)]});
subplot(2,2,2)
plot(D_optim);
xlabel('shell number')
title({'det(F^{-1}) - minimize for D-optimality',['min: ' num2str(D_optim_val)]});
subplot(2,2,3)
plot(E_optim);
xlabel('shell number')
title({'min eigenvalue of F - maximize for E-optimality',['max: ' num2str(E_optim_val)]});
subplot(2,2,4)
plot(T_optim);
xlabel('shell number')
title({'trace(F) - maximiaze for T-optimality',['max: ' num2str(T_optim_val)]});

%% plot eigenvalues progression

figure;
subplot(3,1,1)
plot(eigenvalues(:,1));
[~,optim_val] = max(eigenvalues(:,1));
xlabel('shell number')
title({'F first eigenvalue',['max: ' num2str(optim_val)]})
subplot(3,1,2)
plot(eigenvalues(:,2));
[~,optim_val] = max(eigenvalues(:,2));
xlabel('shell number')
title({'F second eigenvalue',['max: ' num2str(optim_val)]})
subplot(3,1,3)
plot(eigenvalues(:,3));
[~,optim_val] = max(eigenvalues(:,3));
xlabel('shell number')
title({'F third eigenvalue',['max: ' num2str(optim_val)]})

%% Fisher matrix computation

function F = ComputeFisherMatrix(params,bvals,qhat)
DerS0 = DerBallAndStickByS0(params,bvals,qhat);
DerDiff = DerBallAndStickByDiff(params,bvals,qhat);
DerF = DerBallAndStickByF(params,bvals,qhat);

Ders = [DerS0; DerDiff; DerF];

sigma_k = 0.04;

F = zeros(3);

for i=1:3
    for j=1:3
        F(i,j)= (1/sigma_k^2)*sum(Ders(i,:).*Ders(j,:));
    end
end
end


%% derivatives computations functions

function DerS0 = DerBallAndStickByS0(x,bvals,qhat)
% S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
DerS0 = (f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));
end

function DerDiff = DerBallAndStickByDiff(x,bvals,qhat)
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
bqn = bvals.*(fibdotgrad.^2);
DerDiff = S0*(f*(-bqn).*exp(-diff*bqn) + (1-f)*(-bvals).*exp(-bvals*diff));
end

function DerF = DerBallAndStickByF(x,bvals,qhat)
S0 = x(1);
diff = x(2);
% f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
DerF = S0*(exp(-bvals*diff.*(fibdotgrad.^2)) + -1*exp(-bvals*diff));
end




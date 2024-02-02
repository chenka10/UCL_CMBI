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

params = [1.009865516767256   0.001432055891641   0.574928310559141  -1.544730918822476   6.200322564900709];

F = ComputeFisherMatrix(params,bvals,qhat);

%% Find optimal shell
params = [1.009865516767256   0.001432055891641   0.574928310559141  -1.544730918822476   6.200322564900709];

A_optim = zeros(size(shells_idencies,1),1);
D_optim = zeros(size(shells_idencies,1),1);
E_optim = zeros(size(shells_idencies,1),1);
T_optim = zeros(size(shells_idencies,1),1);


for i=1:size(shells_idencies,1)
    curr_shell_idencies = shells_idencies(i,:);
    curr_bvals = bvals(curr_shell_idencies);
    curr_qhat = qhat(curr_shell_idencies);
    F = ComputeFisherMatrix(params,curr_bvals,curr_qhat);
    F_inv = pinv(F);
    A_optim(i)=trace(F_inv);
    D_optim(i) = det(F_inv);
    T_optim(i) = trace(F);

    eigenvalues = eig(F);
    min_eigenvalue = min(eigenvalues);
    E_optim(i) = min_eigenvalue;

end


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




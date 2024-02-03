function [Avox,qhat,TE,bvals] = LoadISBI2015Data(voxel_index)
% Load the diffusion signal
fid = fopen('isbi2015_data_normalised.txt', 'r', 'b');
fgetl(fid); % Read in the header
D = fscanf(fid, '%f', [6, inf])'; % Read in the data
fclose(fid);
% Select the first of the 6 voxels
Avox = D(:,voxel_index);
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
end


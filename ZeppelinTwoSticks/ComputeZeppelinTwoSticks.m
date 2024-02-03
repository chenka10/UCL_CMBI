function S = ComputeZeppelinTwoSticks(x, bvals, qhat)
% Extract the parameters
S0 = x(1);
lam1 = x(2);
lam2 = x(3);
f1 = x(4);
f2 = x(5);
theta_1 = x(6);
phi_1 = x(7);
theta_2 = x(8);
phi_2 = x(9);
theta_3 = x(10);
phi_3 = x(11);

% Synthesize the signals according to the model
fibdir_1 = [cos(phi_1)*sin(theta_1) sin(phi_1)*sin(theta_1) cos(theta_1)];
fibdotgrad_1 = sum(qhat.*repmat(fibdir_1, [length(qhat) 1])');
fibdir_2 = [cos(phi_2)*sin(theta_2) sin(phi_2)*sin(theta_2) cos(theta_2)];
fibdotgrad_2 = sum(qhat.*repmat(fibdir_2, [length(qhat) 1])');
fibdir_3 = [cos(phi_3)*sin(theta_3) sin(phi_3)*sin(theta_3) cos(theta_3)];
fibdotgrad_3 = sum(qhat.*repmat(fibdir_3, [length(qhat) 1])');
S_Z = exp(-bvals.*(lam2 + (lam1 - lam2)*(fibdotgrad_1.^2)));
S_S_1 = exp(-bvals*lam1.*(fibdotgrad_2.^2));
S_S_2 = exp(-bvals*lam1.*(fibdotgrad_3.^2));
S = S0*(f1*S_Z + f2*S_S_1 + (1-f1-f2)*S_S_2);
end
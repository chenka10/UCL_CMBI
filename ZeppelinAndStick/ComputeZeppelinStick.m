function S = ComputeZeppelinStick(x, bvals, qhat)
% Extract the parameters
S0 = x(1);
lam1 = x(2);
lam2 = x(3);
f = x(4);
theta = x(5);
phi = x(6);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S_E = exp(-bvals.*(lam2 + (lam1 - lam2)*(fibdotgrad.^2)));
S_I = exp(-bvals*lam1.*(fibdotgrad.^2));
S = S0*(f*S_I + (1-f)*S_E);
end
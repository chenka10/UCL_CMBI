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
S = S0*(f*exp(-bvals.*(lam2 + (lam1 - lam2)*(fibdotgrad.^2))) + (1-f)*exp(-bvals*lam1));
end
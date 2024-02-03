function S = ComputeZeppelinStickDot(x, bvals, qhat)
% Extract the parameters
S0 = x(1);
lam1 = x(2);
lam2 = x(3);
f1 = x(4);
f2 = x(5);
theta = x(6);
phi = x(7);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S_E = exp(-bvals.*(lam2 + (lam1 - lam2)*(fibdotgrad.^2)));
S_I = exp(-bvals*lam1.*(fibdotgrad.^2));
S = S0*(f1 + f2*S_I + (1-f1-f2)*S_E);
end
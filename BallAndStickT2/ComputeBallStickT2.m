function S = ComputeBallStickT2(x, bvals, qhat,TE)

S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
T2 = 0.07; %[s]

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0*exp(-TE/T2).*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));
% S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));
end
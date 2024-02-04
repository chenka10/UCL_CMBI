function S = ComputeTensorStickDot(x, bvals, qhat)
% Extract the parameters
S0 = x(1);
d_s = x(2);
d_p = x(3);
d_1 = x(4);
d_2 = x(5);
f1 = x(6);
f2 = x(7);
alpha = x(8);
beta = x(9);
gamma = x(10);

% Synthesize the signals according to the model
fibdir_p = [cos(beta)*cos(gamma) cos(beta)*sin(gamma) -sin(beta)];
fibdir_1 = [sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma) sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma) sin(alpha)*cos(beta)];
fibdir_2 = [cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma) cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma) cos(alpha)*cos(beta)];
fibdotgrad_p = sum(qhat.*repmat(fibdir_p, [length(qhat) 1])');
fibdotgrad_1 = sum(qhat.*repmat(fibdir_1, [length(qhat) 1])');
fibdotgrad_2 = sum(qhat.*repmat(fibdir_2, [length(qhat) 1])');
S_E = exp(-bvals.*(d_p*(fibdotgrad_p.^2) + d_1*(fibdotgrad_1.^2) + d_2*(fibdotgrad_2.^2)));
S_I = exp(-bvals*d_s.*(fibdotgrad_p.^2));
S = S0*(f1 + f2*S_I + (1-f1-f2)*S_E);
end
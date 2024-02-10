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
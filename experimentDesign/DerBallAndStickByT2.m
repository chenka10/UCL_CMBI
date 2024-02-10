function DerT2 = DerBallAndStickByT2(x,bvals,qhat,TE)
M0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
T2 = x(6);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
DerT2 = M0*(TE/(T2.^2)).*exp(-TE/T2).*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));
end
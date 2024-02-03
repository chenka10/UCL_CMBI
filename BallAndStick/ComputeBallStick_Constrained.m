function S = ComputeBallStick_Constrained(x, bvals, qhat)
% Extract the parameters
[S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(x);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));
end
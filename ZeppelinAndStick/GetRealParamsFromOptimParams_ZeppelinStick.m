function [x_res] = GetRealParamsFromOptimParams_ZeppelinStick(x)
% [S0,lam1,lam2,f,theta,phi]
x_res = zeros(1,6);

% Extract the parameters
x_res(1) = x(1).^2; % S0
x_res(3) = x(3).^2; % lam2
x_res(2) = x_res(3) + x(2).^2; % lam1
x_res(4) = 1/(1+exp(-x(4))); % f
x_res(5) = x(5); % theta
x_res(6) = x(6); % phi
end

% 
% S0 = x(1);
% lam1 = x(2);
% lam2 = x(3);
% f = x(4);
% theta = x(5);
% phi = x(6);


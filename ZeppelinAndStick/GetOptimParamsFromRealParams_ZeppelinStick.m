function x_res = GetOptimParamsFromRealParams_ZeppelinStick(x)

x_res = zeros(size(x));

x_res(1) = sqrt(x(1));
x_res(3) = sqrt(x(3));
x_res(2) = sqrt(x(2)-x(3));
x_res(4) = log(x(4)/(1-x(4)));
x_res(5) = x(5);
x_res(6) = x(6);
end

% x_res(1) = x(1).^2; % S0
% x_res(3) = x(3).^2; % lam2
% x_res(2) = x_res(2) + x(2).^2; % lam1
% x_res(4) = 1/(1+exp(-x(4))); % f
% x_res(5) = x(5); % theta
% x_res(6) = x(6); % phi

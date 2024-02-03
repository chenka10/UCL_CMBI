function [x_res] = GetRealParamsFromOptimParams_ZeppelinTwoSticks(x)
% [S0,lam1,lam2,f,theta,phi]
x_res = zeros(1,7);

% Extract the parameters
x_res(1) = x(1).^2; % S0>0
x_res(3) = x(3).^2; % lam2>0
x_res(2) = x_res(3) + x(2).^2; % lam1>lam2
x_res(4) = Sigmoid(x(4)); % 0<f1<1
x_res(5) = Sigmoid(x(5))*(1-x_res(4)); % 0<f2<(1-f1)
x_res(6) = x(6); % theta
x_res(7) = x(7); % phi
x_res(8) = x(8); % theta
x_res(9) = x(9); % phi
x_res(10) = x(10); % theta
x_res(11) = x(11); % phi
end


% S0 = x(1);
% lam1 = x(2);
% lam2 = x(3);
% f1 = x(4);
% f2 = x(5);
% theta = x(6);
% phi = x(7);


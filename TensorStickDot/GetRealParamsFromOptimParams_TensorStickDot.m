function [x_res] = GetRealParamsFromOptimParams_TensorStickDot(x)
% [S0,lam1,lam2,f,theta,phi]
x_res = zeros(1,9);

% Extract the parameters
x_res(1) = x(1).^2; % S0>0
x_res(2) = x(2).^2; % d_s>0
x_res(3) = x(3).^2; % d_p>0
x_res(4) = Sigmoid(x(4))*x_res(3); % 0<d_1<d_p
x_res(5) = Sigmoid(x(5))*x_res(4); % 0<d_2<d_1
x_res(6) = Sigmoid(x(6)); % 0<f1<1
x_res(7) = Sigmoid(x(7))*(1-x_res(6)); % 0<f2<(1-f1)
x_res(8) = x(8); % alpha
x_res(9) = x(9); % beta
x_res(10) = x(10); % gamma
end

% S0 = x(1);
% d_p = x(2);
% d_1 = x(3);
% d_2 = x(4);
% f1 = x(5);
% f2 = x(6);
% alpha = x(7);
% beta = x(8);
% gamma = x(9);


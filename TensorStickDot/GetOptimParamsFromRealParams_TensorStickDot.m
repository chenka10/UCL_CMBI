function x_res = GetOptimParamsFromRealParams_TensorStickDot(x)

x_res = zeros(size(x));

x_res(1) = sqrt(x(1));
x_res(2) = sqrt(x(2));
x_res(3) = sqrt(x(3));
x_res(4) = SigmoidInv(x(4)/x(3));
x_res(5) = SigmoidInv(x(5)/x(4));
x_res(6) = SigmoidInv(x(6));
x_res(7) = SigmoidInv(x(7)/(1-x(6)));
x_res(8) = x(8);
x_res(9) = x(9);
x_res(10) = x(10);
end

% Extract the parameters
% x_res(1) = x(1).^2; % S0>0
% x_res(2) = x(2).^2; % d_p>0
% x_res(3) = Sigmoid(x(3))*x_res(2); % 0<d_1<d_p
% x_res(4) = Sigmoid(x(4))*x_res(3); % 0<d_2<d_1
% x_res(5) = Sigmoid(x(5)); % 0<f1<1
% x_res(6) = Sigmoid(x(6))*(1-x_res(5)); % 0<f2<(1-f1)
% x_res(7) = x(7); % alpha
% x_res(8) = x(8); % beta
% x_res(9) = x(9); % gamma

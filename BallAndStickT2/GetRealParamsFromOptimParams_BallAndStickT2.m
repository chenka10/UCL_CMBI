function x_res = GetRealParamsFromOptimParams_BallAndStickT2(x)

x_res = zeros(6,1);

% Extract the parameters
x_res(1) = x(1).^2;
x_res(2) = x(2).^2;
x_res(3) = 1/(1+exp(-x(3)));
x_res(4) = x(4);
x_res(5) = x(5);
x_res(6) = x(6).^2;
end


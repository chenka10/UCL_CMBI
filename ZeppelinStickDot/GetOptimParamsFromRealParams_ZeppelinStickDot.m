function x_res = GetOptimParamsFromRealParams_ZeppelinStickDot(x)

x_res = zeros(size(x));

x_res(1) = sqrt(x(1));
x_res(3) = sqrt(x(3));
x_res(2) = sqrt(x(2)-x(3));
x_res(4) = SigmoidInv(x(4));
x_res(5) = SigmoidInv(x(5)/(1-x(4)));
x_res(6) = x(6);
x_res(7) = x(7);
end

% x_res(1) = x(1).^2; % S0>0
% x_res(3) = x(3).^2; % lam2>0
% x_res(2) = x_res(3) + x(2).^2; % lam1>lam2
% x_res(4) = Sigmoid(x(4)); % 0<f1<1
% x_res(5) = Sigmoid(x(5))*(1-x_res(4)); % 0<f2<(1-f1)
% x_res(6) = x(6); % theta
% x_res(7) = x(7); % phi

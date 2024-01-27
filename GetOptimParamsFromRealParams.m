function x_res = GetOptimParamsFromRealParams(x)

x_res = zeros(size(x));

x_res(1) = sqrt(x(1));
x_res(2) = sqrt(x(2));
x_res(3) = log(x(3)/(1-x(3)));
% x_res(4) = tan(x(4)-pi/2);
x_res(4) = x(4);
x_res(5) = x(5);
end


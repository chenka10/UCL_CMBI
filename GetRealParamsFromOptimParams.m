function [S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(x)
% Extract the parameters
S0 = x(1).^2;
diff = x(2).^2;
f = 1/(1+exp(-x(3)));
% theta = tanh(x(4))+pi/2;
theta = x(4);
phi = x(5);
end


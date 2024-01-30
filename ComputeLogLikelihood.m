function [L] = ComputeLogLikelihood(sigma,Data,ModelSignal)
components = -0.5*log(2*pi*(sigma^2))-((Data - ModelSignal).^2)/(2*(sigma^2));
% (1/sqrt(2*pi*(sigma^2)))*exp((-(Data - ModelSignal).^2)/(2*(sigma^2)));
L = sum(components);
end


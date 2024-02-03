function AIC = ComputeAIC(N,K,SSD)
AIC = 2*N + K*log((1/K)*SSD);
end
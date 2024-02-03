
function BIC = ComputeBIC(N,K,SSD)
BIC = N*log(K) + K*log((1/K)*SSD);
end
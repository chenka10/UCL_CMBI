clc;
K=3612;
Ns = [8 6 7 6];

for i=1:numel(Ns)
    N = Ns(i);
    SSD = sum((Avox-models_res{i}').^2);
    AIC = ComputeAIC(N,K,SSD);
    BIC = ComputeBIC(N,K,SSD);

    fprintf('%s; AIC: %d BIC: %d\n',names{i},AIC,BIC);
end



function F = ComputeFisherMatrixWithT2(params,bvals,qhat,TE)
T2 = params(6);
T2_exp_component = exp(-TE/T2);

params = [params T2];

DerM0 = T2_exp_component.*DerBallAndStickByS0(params,bvals,qhat);
DerDiff = T2_exp_component.*DerBallAndStickByDiff(params,bvals,qhat);
DerF = T2_exp_component.*DerBallAndStickByF(params,bvals,qhat);
DerT2 = DerBallAndStickByT2(params,bvals,qhat,TE);


Ders = [DerM0; DerDiff; DerF; DerT2];

sigma_k = 0.04;

F = zeros(4);

for i=1:4
    for j=1:4
        F(i,j)= (1/sigma_k^2)*sum(Ders(i,:).*Ders(j,:));
    end
end
end
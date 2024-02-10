function F = ComputeFisherMatrix(params,bvals,qhat,TE)
DerS0 = DerBallAndStickByS0(params,bvals,qhat);
DerDiff = DerBallAndStickByDiff(params,bvals,qhat);
DerF = DerBallAndStickByF(params,bvals,qhat);

Ders = [DerS0; DerDiff; DerF];

sigma_k = 0.04;

F = zeros(3);

for i=1:3
    for j=1:3
        F(i,j)= (1/sigma_k^2)*sum(Ders(i,:).*Ders(j,:));
    end
end
end
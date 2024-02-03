function sumRes = SSD_DTI(x, Avox, bvals, qhat)

S = ComputeDti(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
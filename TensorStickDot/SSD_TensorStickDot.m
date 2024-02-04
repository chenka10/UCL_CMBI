function sumRes = SSD_TensorStickDot(x, Avox, bvals, qhat)

x = GetRealParamsFromOptimParams_TensorStickDot(x);

S = ComputeTensorStickDot(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
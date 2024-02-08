function sumRes = BallStickT2_SSD(x, Avox, bvals, qhat, TE)

% Extract the parameters
x = GetRealParamsFromOptimParams_BallAndStickT2(x);


S = ComputeBallStickT2(x, bvals, qhat, TE);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
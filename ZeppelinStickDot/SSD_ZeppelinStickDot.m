function sumRes = SSD_ZeppelinStickDot(x, Avox, bvals, qhat)

x = GetRealParamsFromOptimParams_ZeppelinStickDot(x);

S = ComputeZeppelinStickDot(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
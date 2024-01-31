function sumRes = SSD_ZeppelinAndStick(x, Avox, bvals, qhat)

x = GetRealParamsFromOptimParams_ZeppelinStick(x);

S = ComputeZeppelinStick(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
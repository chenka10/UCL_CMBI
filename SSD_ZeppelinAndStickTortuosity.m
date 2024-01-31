function sumRes = SSD_ZeppelinAndStickTortuosity(x, Avox, bvals, qhat)

[S0,diff,f,theta,phi] = GetRealParamsFromOptimParams(x);
x = [S0,diff,f,theta,phi];


S = ComputeZeppelinStickTortuosity(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
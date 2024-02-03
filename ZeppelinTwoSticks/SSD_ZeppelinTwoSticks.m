function sumRes = SSD_ZeppelinTwoSticks(x, Avox, bvals, qhat)

x = GetRealParamsFromOptimParams_ZeppelinTwoSticks(x);

S = ComputeZeppelinTwoSticks(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
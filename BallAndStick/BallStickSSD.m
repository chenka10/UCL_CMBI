function sumRes = BallStickSSD(x, Avox, bvals, qhat)
S = ComputeBallStick(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
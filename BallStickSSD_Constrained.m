function sumRes = BallStickSSD_Constrained(x, Avox, bvals, qhat)

S = ComputeBallStick_Constrained(x, bvals, qhat);

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);

end
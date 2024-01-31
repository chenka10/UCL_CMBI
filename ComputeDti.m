function S = ComputeDti(x, bvals, qhat)

D = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];
S0 = x(1);

qDq = sum(qhat.*(D*qhat));

S = S0*(exp(-bvals.*qDq));
end
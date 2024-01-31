function D = GetDiffusionMatrixForVoxel(Avox,Y_pinv)

x = Y_pinv*log(Avox);

D = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];
end
function [starting_values,fitted_params,resnorms] = RandomBallStickFitting(startx_orig,noise_range,Avox,qhat,bvals,N)

starting_values = zeros(N,5);
fitted_params = zeros(N,5);
resnorms = zeros(N,1);

for i=1:N

startx = startx_orig;

noise = randn(1,5).*noise_range;
startx = startx + noise;
startx(1) = abs(startx(1));
startx(2) = abs(startx(2));
startx(3) = max(min(startx(3),1),0);
startx(4) = max(min(startx(4),pi),0);

starting_values(i,:) = startx;

startx(1) = sqrt(startx(1));
startx(2) = sqrt(startx(2));
startx(3) = log(startx(3)/(1-startx(3)));
startx(4) = tan(startx(4)-pi/2);

h=optimset('MaxFunEvals',20000,...
 'Algorithm','quasi-newton',...
 'MaxIter',200,...
 'TolX',1e-10,...
 'TolFun',1e-10, ...
 'Display','none');

[parameter_hat,RESNORM,~,~]=fminunc('BallStickSSD_Constrained',startx,h,Avox,bvals,qhat);

fitted_params(i,:) = parameter_hat;
resnorms(i) = RESNORM;

end
end


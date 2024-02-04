function [starting_values,fitted_params, resnorms, hessians] = RandomBallStickFitting(startx_orig,noise_range,Avox,qhat,bvals,N)

starting_values = zeros(N,5);
fitted_params = zeros(N,5);
resnorms = zeros(N,1);
hessians = cell(N,1);

for i=1:N
    try
        startx = startx_orig;

        % generate randomness
        noise = randn(1,5).*noise_range;
        startx = startx + noise;

        % clip perturbed params
        startx(1) = abs(startx(1));
        startx(2) = abs(startx(2));
        startx(3) = mod(startx(3),1);        

        % save start values
        starting_values(i,:) = startx;
        
        % inverse start params
        startx = GetOptimParamsFromRealParams(startx);

        h=optimset('MaxFunEvals',20000,...
            'Algorithm','quasi-newton',...
            'MaxIter',50,...
            'TolX',1e-10,...
            'TolFun',1e-10, ...
            'Display','none');

        [parameter_hat,RESNORM,~,~,~,hessian]=fminunc('BallStickSSD_Constrained',startx,h,Avox,bvals,qhat);

        fitted_params(i,:) = parameter_hat;
        resnorms(i) = RESNORM;
        hessians{i} = hessian;
    catch
        resnorms(i) = inf;
    end
end
end


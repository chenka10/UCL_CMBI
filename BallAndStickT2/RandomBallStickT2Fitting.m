function [starting_values,fitted_params, resnorms, hessians] = RandomBallStickT2Fitting(startx_orig,noise_range,Avox,qhat,bvals,TE,N)

starting_values = zeros(N,6);
fitted_params = zeros(N,6);
resnorms = zeros(N,1);
hessians = cell(N,1);

for i=1:N
    try
        startx = startx_orig;

        % generate randomness
        noise = randn(1,6).*noise_range;
        startx = startx + noise;

        % clip perturbed params
        startx(1) = abs(startx(1));
        startx(2) = abs(startx(2));
        startx(3) = mod(startx(3),1);    
        startx(6) = mod(startx(6),1);

        % save start values
        starting_values(i,:) = startx;
        
        % inverse start params
        startx = GetOptimParamsFromRealParams_BallAndStickT2(startx);

        h=optimset('MaxFunEvals',20000,...
            'Algorithm','quasi-newton',...
            'MaxIter',500,...
            'TolX',1e-10,...
            'TolFun',1e-10, ...
            'Display','none');

        [parameter_hat,RESNORM,~,~,~,hessian]=fminunc('BallStickT2_SSD',startx,h,Avox,bvals,qhat,TE);

        fitted_params(i,:) = parameter_hat;
        resnorms(i) = RESNORM;
        hessians{i} = hessian;
    catch
        resnorms(i) = inf;
    end
end
end


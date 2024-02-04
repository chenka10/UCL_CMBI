function [starting_values,fitted_params, resnorms] = RandomTensorStickDotFitting(startx_orig,noise_range,Avox,qhat,bvals,N)


starting_values = zeros(N,10);
fitted_params = zeros(N,10);
resnorms = zeros(N,1);

for i=1:N
    try
        startx = startx_orig;

        % generate randomness
        noise = randn(1,10).*noise_range;
        startx = startx + noise;

        % clip perturbed params
        startx(1) = abs(startx(1));        
        startx(2) = abs(startx(2));
        startx(3) = abs(startx(3));
        startx(4) = min(startx(3)*0.9999,abs(startx(4)));
        startx(5) = min(startx(4)*0.9999,abs(startx(5)));
        startx(6) = mod(startx(6),1);    
        startx(7) = mod(startx(7),1)*(1-startx(6));

        % save start values
        starting_values(i,:) = startx;
        
        % inverse start params
        startx = GetOptimParamsFromRealParams_TensorStickDot(startx);

        h=optimset('MaxFunEvals',20000,...
            'Algorithm','quasi-newton',...
            'MaxIter',200,...
            'TolX',1e-10,...
            'TolFun',1e-10, ...
            'Display','none');

        [parameter_hat,RESNORM,~,~]=fminunc('SSD_TensorStickDot',startx,h,Avox,bvals,qhat);

        fitted_params(i,:) = parameter_hat;
        resnorms(i) = RESNORM;
    catch
        resnorms(i) = inf;
    end
end
end


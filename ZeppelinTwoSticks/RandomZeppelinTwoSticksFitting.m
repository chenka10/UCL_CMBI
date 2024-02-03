function [starting_values,fitted_params, resnorms] = RandomZeppelinTwoSticksFitting(startx_orig,noise_range,Avox,qhat,bvals,N)


starting_values = zeros(N,11);
fitted_params = zeros(N,11);
resnorms = zeros(N,1);

for i=1:N
    try
        startx = startx_orig;

        % generate randomness
        noise = randn(1,11).*noise_range;
        startx = startx + noise;

        % clip perturbed params
        startx(1) = abs(startx(1));        
        startx(3) = abs(startx(3));
        startx(2) = max(startx(2),startx(3));
        startx(4) = mod(startx(4),1);    
        startx(5) = mod(startx(5),1)*(1-startx(4));

        % save start values
        starting_values(i,:) = startx;
        
        % inverse start params
        startx = GetOptimParamsFromRealParams_ZeppelinTwoSticks(startx);

        h=optimset('MaxFunEvals',20000,...
            'Algorithm','quasi-newton',...
            'MaxIter',200,...
            'TolX',1e-10,...
            'TolFun',1e-10, ...
            'Display','none');

        [parameter_hat,RESNORM,~,~]=fminunc('SSD_ZeppelinTwoSticks',startx,h,Avox,bvals,qhat);

        fitted_params(i,:) = parameter_hat;
        resnorms(i) = RESNORM;
    catch
        resnorms(i) = inf;
    end
end
end


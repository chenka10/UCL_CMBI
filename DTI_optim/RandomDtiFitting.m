function [starting_values,fitted_params, resnorms] = RandomDtiFitting(startx_orig,noise_range,Avox,qhat,bvals,N)

starting_values = zeros(N,7);
fitted_params = zeros(N,7);
resnorms = zeros(N,1);

for i=1:N
    try
        startx = startx_orig;

        % generate randomness
        noise = randn(1,7).*noise_range;
        startx = startx + noise;

        % save start values
        starting_values(i,:) = startx;        

        h=optimset('MaxFunEvals',20000,...
            'Algorithm','quasi-newton',...
            'MaxIter',200,...
            'TolX',1e-10,...
            'TolFun',1e-10, ...
            'Display','none');

        [parameter_hat,RESNORM,~,~]=fminunc('SSD_DTI',startx,h,Avox,bvals,qhat);

        fitted_params(i,:) = parameter_hat;
        resnorms(i) = RESNORM;
    catch
        resnorms(i) = inf;
    end
end
end


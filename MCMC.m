
function [results,acceptance_count] = MCMC(Avox,qhat,bvals,Y_pinv,N,stabilization_iterations)

x = Y_pinv*log(Avox);

D_DTI = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];
S0_DTI = exp(x(1));

R = D_DTI/trace(D_DTI);
FA_DTI = sqrt(0.5*(3-1/trace(R^2)));

N = N+stabilization_iterations;

% we assume standard deviation of 200
sigma = 200;

results = zeros(N,6);

% acceptance_log
acceptance_count = 0;

% setup random noise range to fit parameter values
S0_range = 30;
d_range = 0.0001;
f_range = 0.003;
theta_range = pi/100;
phi_range = pi/100;
noise_range = [S0_range, d_range, f_range, theta_range, phi_range];

% setup starting params
current_params = [S0_DTI trace(D_DTI)/3 FA_DTI 0 0];

current_model_signals = ComputeBallStick(current_params,bvals,qhat)';
current_log_likelihood = ComputeLogLikelihood(sigma,Avox,current_model_signals);

for i=1:N

    % perturbe params
    noise = randn(1,5).*noise_range;
    params = current_params + noise;

    % limit perturbed params
    params(1) = abs(params(1));
    params(2) = abs(params(2));
    params(3) = mod(params(3),1);

    current_model_signals = ComputeBallStick(params,bvals,qhat)';

    resnorm = sum((Avox - current_model_signals).^2);
    log_likelihood = ComputeLogLikelihood(sigma,Avox,current_model_signals);

    if (exp(log_likelihood-current_log_likelihood) > rand())
        current_params = params;
        current_log_likelihood = log_likelihood;
        if (i>stabilization_iterations)
            acceptance_count = acceptance_count + 1;
        end
    end

    results(i,:) = [current_params,resnorm];

    if mod(i,10000)==0
        disp([num2str(i) '/' num2str(N)])
    end

end
end
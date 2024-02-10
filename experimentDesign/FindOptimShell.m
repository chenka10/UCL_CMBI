function [A_optim,D_optim,E_optim,T_optim,eigenvalues] = FindOptimShell(shells_idencies, params,params_norm,computeFisherMatrixFunc,bvals,qhat,TE)

A_optim = zeros(size(shells_idencies,1),1);
D_optim = zeros(size(shells_idencies,1),1);
E_optim = zeros(size(shells_idencies,1),1);
T_optim = zeros(size(shells_idencies,1),1);

for i=1:size(shells_idencies,1)
    curr_shell_idencies = shells_idencies(i,:);
    curr_bvals = bvals(curr_shell_idencies);
    curr_qhat = qhat(curr_shell_idencies);
    curr_TE = TE(curr_shell_idencies);
    fprintf("shell %d has TE %d\n",i,curr_TE(1))
    F = computeFisherMatrixFunc(params,curr_bvals,curr_qhat,curr_TE).*params_norm;
    F_inv = pinv(F);
    A_optim(i)= trace(F_inv);
    D_optim(i) = det(F_inv);

    curr_eigenvalues = eig(F);
    min_eigenvalue = min(curr_eigenvalues);
    E_optim(i) = min_eigenvalue;

    if i==1
        eigenvalues = zeros(size(shells_idencies,1),numel(curr_eigenvalues));
    end

    eigenvalues(i,:) = curr_eigenvalues;

    T_optim(i) = F(1,1);
end
end


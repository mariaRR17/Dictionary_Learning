function [D_new, X_iter, sp, res,mu] = MOD_iterative_th( D_new, X_iter, Y, mu, lambda, numiter_Sparse_Coding,  numiter_dictionary)

 for j=1:numiter_dictionary
    
        %Sparse coding step: ISTA with gradient step mu 
        X_iter = iterative_soft_thresholding(D_new,X_iter,Y,mu,lambda,numiter_Sparse_Coding);
    
        
        % MOD algorithm
        D_new = Y * pinv(X_iter);
        D_new = normalize_columns(D_new);
        mu = 1.9/norm(D_new)^2;
        
        sp(j)=sum(sum(abs(X_iter)>1e-5)); % to monitor sparsity decrease
        res(j) = norm(Y-D_new*X_iter,'fro');
 end  



end


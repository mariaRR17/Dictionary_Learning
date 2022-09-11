function [D_new, X_iter, sp, res] = MOD_BPDN( D_new, X_iter, Y, lambda,  numiter_dictionary,l)

 for j=1:numiter_dictionary
    
        %Sparse coding step: basis pursuit denoinsing with parameter lambda
        for i=1:l
        X_iter(:,i) = basis_pursuit_DN(D_new,Y(:,i), lambda);  
         
        end
        
        % MOD algorithm
        D_new = Y * pinv(X_iter);
        D_new = normalize_columns(D_new);
        
        sp(j)=sum(sum(abs(X_iter)>1e-5)); % to monitor sparsity decrease
        res(j) = norm(Y-D_new*X_iter,'fro');
 end  



end


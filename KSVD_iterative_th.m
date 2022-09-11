function [D_new, X_iter, sp, res,mu] = KSVD_iterative_th( D_new, X_iter, Y, mu, lambda, numiter_Sparse_Coding,  numiter_dictionary,m, k )
     
     for j=1:numiter_dictionary
    
        %Sparse coding step: ISTA with gradient step mu 
        X_iter = iterative_soft_thresholding(D_new,X_iter,Y,mu,lambda,numiter_Sparse_Coding);
    
        
        % K-SVD algorithm
    
        for a = 1:m
            y_bar = [];
            x_bar = [];
            omega = [];
            for i = 1:k
                if abs(X_iter(a,i)) > 1e-5
                    %signals that have used atom j
                    y_bar = [y_bar,Y(:,i)];
                    x_bar = [x_bar, X_iter(:,i)];
                    omega = [omega,i];
    
                end 
    
            end
            if ~isempty(y_bar)
                %contributions of all atoms to the signal y, except atom j
                D_bar = D_new;
                D_bar(:,a) = [];
                x_bar(a,:) = [];
                %Remove all the contributions of all atoms, but the chose one:j
                e = y_bar- D_bar*x_bar;
                %Apply SVD to e: e = U*S*V'
                [U,S,V] = svd(e);
                D_new(:,a) = U(:,1);
                X_iter(a,omega) = S(1,1)*V(:,1)';
    
            end
    
        end
        %D_new = normalize_columns(D_new);
        mu = 1.9/norm(D_new)^2;
        
        sp(j)=sum(sum(abs(X_iter)>1e-5)); % to monitor sparsity decrease
        res(j) = norm(Y-D_new*X_iter,'fro');
     end    
end  






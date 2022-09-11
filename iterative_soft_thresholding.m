

function [X_iter] = iterative_soft_thresholding(D_new,X_iter,Y,mu,lambda,numiter_Sparse_Coding)


for i=1:numiter_Sparse_Coding
    X_iter = X_iter + mu*D_new'*(Y-D_new*X_iter);  % gradient step
    % soft thresholding
    signo = sign(X_iter);
    X_iter = abs(X_iter) - lambda*mu;
    X_iter = (X_iter + abs(X_iter))/2;
    X_iter = signo.*X_iter;
    
    %to visualize the inner iterations to calibrate lambda and numiter_Sparse_Coding
    %plot(X_iter(:,1),'-o')
    %pause(0.1)
    %i
end

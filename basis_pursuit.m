%Basis pursuit algorithm

function [x] = basis_pursuit(D,y)
options = optimoptions('linprog','Display','none');
%Normalize D 
[m,n] = size(D);
W = eye(n);
for i = 1:n
    column_norm = norm(D(:,i));
    W(i,i) = 1/column_norm;
end

D_bar = D*W;

%Solve linear progamming problem
%Objective function
f = [ zeros(1,n), ones(1,n)];

%Inequality constraints
I = eye(n);
A = [I, -I; -I, -I];
b = zeros(1,2*n);

%Equality constraints
Aeq = [D_bar,zeros(m,n)];
beq = y;
x_u_bar = linprog(f,A,b,Aeq,beq,[],[],options);
x_bar = x_u_bar(1:n);

%Denormalized x
x = W*x_bar;

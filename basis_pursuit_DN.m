%Basis pursuit algorithm

function [x] = basis_pursuit_DN(D,y, lambda)
options = optimoptions('quadprog','Display','none');
%Normalize D 
[m,n] = size(D);
W = eye(n);
for i = 1:n
    column_norm = norm(D(:,i));
    W(i,i) = 1/column_norm;
end

D_bar = D*W;

%Solve quadratic progamming problem
%Objective function
Q= [D_bar,-D_bar];
H = 2*(Q'*Q);

c =lambda*[ ones(n,1); ones(n,1)];

f = -2*Q'*y + c;
%Inequality constraints
I = eye(n);
A = [];
b = [];

%Equality constraints
Aeq = [];
beq = [];
lb = zeros(2*n,1);
x_u_bar = quadprog(H,f,A,b,Aeq,beq,lb,[],[],options);
u_bar = x_u_bar(1:n);
v_bar = x_u_bar(n+1:2*n);

%Denormalized x
x_bar = u_bar-v_bar;
x = W*x_bar;

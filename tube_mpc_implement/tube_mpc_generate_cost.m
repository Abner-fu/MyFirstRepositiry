function oc_cost = tube_mpc_generate_cost(cost, N)

n=size(cost.Q,1);
m=size(cost.R,1);

H=zeros((N+1)*n+N*m,(N+1)*n+N*m);
H(1:N*n,1:N*n)=kron(eye(N),cost.Q);
H(N*n+1:(N+1)*n,N*n+1:(N+1)*n)=cost.P;
H((N+1)*n+1:(N+1)*n+N*m,(N+1)*n+1:(N+1)*n+N*m)=kron(eye(N),cost.R);
% take care of QP format used by quadprog
oc_cost.H=2*H; 
oc_cost.f=zeros((N+1)*n+N*m,1);

end
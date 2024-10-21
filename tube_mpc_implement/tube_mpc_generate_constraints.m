% 注意，这里的控制变量是（z_0, z_1, z_2, ..., z_N-1, z_N, v_0, v_1, v_2, ..., v_N-1）
% 不等式约束 z_k 属于 Z
%           v_k 属于 V
%           z_N 属于 Zf
%           x-z_0 属于 S

% 等式约束  z_k+1 = A*z_k + B*v_k
function oc_constraints = tube_mpc_generate_constraints(nominal_system, Z, V, Zf, S, N, x)
    % 生成不等式约束 Ain 和 bin
    n_S = size(S.A);
    n_Z = size(Z.A);
    n_V = size(V.A);
    n_Zf = size(Zf.A);
    oc_constraints.Ain = zeros( n_Z(1)*N + n_Zf(1) + n_V(1)*N + n_S(1), n_Z(2)*N + n_Zf(2) + n_V(2)*N);
    oc_constraints.Ain(1:n_Z(1)*N, 1:n_Z(2)*N) = kron(eye(N), Z.A);
    oc_constraints.Ain(n_Z(1)*N+1:n_Z(1)*N+n_Zf(1), n_Z(2)*N+1:n_Z(2)*N+n_Zf(2)) = Zf.A;
    oc_constraints.Ain(n_Z(1)*N+n_Zf(1)+1:n_Z(1)*N+n_Zf(1)+n_V(1)*N, n_Z(2)*N+n_Zf(2)+1:n_Z(2)*N+n_Zf(2)+n_V(2)*N) = kron(eye(N), V.A);
    oc_constraints.Ain(n_Z(1)*N+n_Zf(1)+n_V(1)*N+1:n_Z(1)*N+n_Zf(1)+n_V(1)*N+n_S(1), 1:n_S(2)) = - S.A;

    oc_constraints.bin = zeros(n_Z(1)*N + n_Zf(1) + n_V(1)*N + n_S(1), 1);
    oc_constraints.bin(1:n_Z(1)*N, 1) = kron(ones(N, 1), Z.b);
    oc_constraints.bin(n_Z(1)*N+1:n_Z(1)*N+n_Zf(1), 1) = Zf.b;
    oc_constraints.bin(n_Z(1)*N+n_Zf(1)+1:n_Z(1)*N+n_Zf(1)+n_V(1)*N, 1) = kron(ones(N, 1), V.b);
    oc_constraints.bin(n_Z(1)*N+n_Zf(1)+n_V(1)*N+1:n_Z(1)*N+n_Zf(1)+n_V(1)*N+n_S(1), 1) = S.b - (S.A * x);
    
    n = size(nominal_system.A, 2);
    m = size(nominal_system.B, 2);
    oc_constraints.Aeq = zeros(n*N, (N+1)*n+N*m);
    for i = 0:N-1
        oc_constraints.Aeq(1+i*n:n+i*n, 1+n*i:n+n*i) = -nominal_system.A;
        oc_constraints.Aeq(1+i*n:n+i*n, n+1+n*i:2*n+n*i) = eye(n);
        oc_constraints.Aeq(1+i*n:n+i*n, (N+1)*n+1+m*i:(N+1)*n+m+m*i) = -nominal_system.B;
    end
    
    oc_constraints.beq = zeros(n*N, 1);
end
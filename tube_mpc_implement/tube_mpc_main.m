% 在我的实现里和论文有些符号不太一样， 论文里的 Z 对应代码里的 S
%                                          Xf 对应代码里的 Zf
% 这里 nominal system 符号使用的是 z+ = A*z + B*v
% 这里 system 符号是 x+ = A*x + B*u + w

% 实现tube mpc 是先离线构造出需要的集合, 然后再在线计算，
% 下面是需要离线构造的集合:
% control feedback matrix: K
% Perturbation Invariant Set(minimal robust positive invariant set): S
% Control Constraint Fields for Nominal Systems:  V = U - K*S
% Domain of the nominal system： Z = X - S
% target set of the nominal system：Zf 

% 当控制系统到达Zf后，使他不会跑出这个集合，采用的控制使用的是 v = Kz , 
% 所以在求 Zf【一般取 nominal system（此时是z+=(A+B*K)z）的最大正不变集】的
% 时候 z 需要满足的约束是 { z: z belongs to Z, K*z belongs to V }

close all;
clear
clc

rng(6)

A = [1, 1; 0, 1];
B = [0.5; 1];

Q = (1/2)*eye(2);
R = (1/2)*0.01;

[K, P] = dlqr(A, B, Q, R);
K = -K;

s_system.A = A+B*K;
s_system.b = [0; 0];

disturbance.E = [1, 0; -1, 0; 0, 1; 0, -1];
disturbance.g = [0.1; 0.1; 0.1; 0.1];

% take alpha in (0, 1), The smaller the better but too small maybe no solution.
alpha = 0.1;

[I, S] = Minimal_robust_positively_invariant_set(s_system, disturbance, alpha);
X = Polyhedron([0, 1], 2);
U = Polyhedron([1; -1], [1; 1]);
% disp(S <= X)
% disp(K*S <= U)
% 
% figure()
% plot(S)

Z = X - S;
V = U - K*S;
% figure()
% plot([X, Z])
% figure()
% plot([U, V])

z_system.A = A+B*K;
z_system.b = [0; 0];

z_constraints.C = [Z.A; V.A * K]; 
z_constraints.d = [Z.b; V.b];

Zf = Mpis(z_system, z_constraints, 100);
% figure();
% plot([Zf+S Zf])

cost_param.Q = Q;
cost_param.P = P;
cost_param.R = R;
N = 9;

nominal_system.A = A;
nominal_system.B = B;
x0 = [-5; -2];

% 方法一，只是解决一次优化问题
% N = 9;
% optimal = tube_mpc_oc(cost_param, nominal_system, Z, V, Zf, S, N, x0);
% v = optimal.v;
% z = optimal.z;
% x(:, 1) = x0;
% for i = 1:N
%     u(:, i) = v(:, i) + K*(x(:, i) - z(:, i));
%     w(:, i) = rand()*0.2 - 0.1;
%     x(:, i+1) = A*x(:, i) + B*u(:, i) + w(:, i);
% end
% 
% % 数据处理，以便绘图
% tube = [];
% for i = 1:N+1
%     tube = [tube; z(:, i) + S];
% end
% figure()
% plot([Zf+S Zf tube]); hold on;
% plot(x(1, :), x(2, :), '--'); hold on;
% plot(z(1, :), z(2, :), '-')

% 方法二，每个时刻都解一次优化问题
% tube mpc online construct

Nsim = 9;
x(:, 1) = x0;
for i = 1:Nsim
    optimal(i) = tube_mpc_oc(cost_param, nominal_system, Z, V, Zf, S, N, x(:, i));
    v(:, i) = optimal(i).v(:, 1);
    z(:, i) = optimal(i).z(:, 1);
    u(:, i) = v(:, i) + K*(x(:, i) - z(:, i));
    w(:, i) = rand() * 0.2 - 0.1;
    % w(:, i) = -0.1;
    x(:, i+1) = A*x(:, i) + B*u(:, i) + w(:, i);
    % N = N-1;
end


% 数据处理，以便绘图
tube = [];
for i = 1:Nsim
    tube = [tube; z(:, i) + S];
end
figure()
plot([Zf+S Zf tube]); hold on;
plot(x(1, 1:Nsim), x(2, 1:Nsim), '--'); hold on;
plot(z(1, :), z(2, :), '-')
xlim([-8, 4]); ylim([-3, 3]);
xlabel("x_1"); ylabel("x_2");

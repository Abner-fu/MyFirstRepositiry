function optimal = tube_mpc_oc(cost_param, nominal_system, Z, V, Zf, S, N, x)
    n = size(nominal_system.A, 2);
    m = size(nominal_system.B, 2);

    oc_cost = tube_mpc_generate_cost(cost_param, N);
    mpc_constraints = tube_mpc_generate_constraints(nominal_system, Z, V, Zf, S, N, x);

    [optimal.solution, optimal.V] = quadprog(oc_cost.H, oc_cost.f, mpc_constraints.Ain, mpc_constraints.bin, mpc_constraints.Aeq, mpc_constraints.beq);
    z(:, 1) = optimal.solution(1:n);
    for i = 1:N
        z(:, i+1) = optimal.solution(i*n+1:(i+1)*n);
        v(:, i) = optimal.solution((N+1)*n+(i-1)*m+1:(N+1)*n+i*m);
    end
    optimal.v = v;
    optimal.z = z;
end
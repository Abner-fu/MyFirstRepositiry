function Z = One_step_backward_reachable_set(system, constraints, target)
M = [constraints.C; target.G*system.A];
L = [constraints.d; target.h - target.G * system.b];

Z = Polyhedron(M, L);
end
function [X_inner, X_outer] = Minimal_robust_positively_invariant_set(system, disturbance, alpha)

i = 1;
W = Polyhedron(disturbance.E, disturbance.g);
X(i) = W;
while not(system.A ^ i * W <= alpha*W)
    i = i + 1;
    X(i) = system.A * X(i-1) + W;
end
x_dim = size(system.A, 2);
X_translation = inv(eye(x_dim) - system.A) * system.b;
X_outer = X_translation + 1/(1-alpha) * X(i);
X_inner = X_translation + X(i);

end
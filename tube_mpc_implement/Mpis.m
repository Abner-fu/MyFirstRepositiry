% 计算最大正不变集（没有控制 u 的情况）
function [X_mpis, X, i] = Mpis(system, constraints, imax)
    i = 1;
    X(i) = Polyhedron(constraints.C, constraints.d);
    cond = true;
    while cond
        i = i + 1;
        current_target.G = X(i-1).A;
        current_target.h = X(i-1).b;
        X(i) = One_step_backward_reachable_set(system, constraints, current_target);
        cond = and(not(X(i-1)==X(i)), i <= imax);
    end

    if i > imax
        disp("cannot find mpis set in " + imax + " steps");
        return;
    end
    X_mpis = X(i);
end
function x = solve_velocity_system(S,b)
%SOLVE_VELOCITY_SYSTEM Solve K*x=b using the stored, unchanged LU factors.
if isempty(b) 
    x = b; 
    return; 
end
F = S.factor;
x = F.Q*(F.U\(F.L\(F.P*b)));
end

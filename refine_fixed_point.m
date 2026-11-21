function [S,info,first] = refine_fixed_point(S,tol,maxiter)
%REFINE_FIXED_POINT Fixed K iteration for actual C2 acceleration matching.
%   K*c^(k+1) = f - N(c^k)
% is evaluated in its equivalent residual form
%   K*delta = J(c^k)/2;  c^(k+1) = c^k + delta.
% J is the RIGHT-minus-LEFT actual covariant acceleration jump. The plus
% sign and the factor 1/2 follow from J=-2*(K*c-f+N).
% The LU factors stored by interlocking_spline are reused, without a new
% Jacobian or a line search. Only the limit is exactly C2. A finite return
% marked converged satisfies the requested actual-jump tolerance.
% Failure is returned explicitly; a failed iterate is not labeled C2.
% FIRST records c^1 when a correction is performed; if c^0 already meets
% the tolerance, FIRST is the unmodified initializer and info.iterations=0.
if nargin < 2 
    tol = 1e-9; 
end
if nargin < 3 
    maxiter = 100; 
end
if ~isscalar(tol) || ~isfinite(tol) || tol <= 0 || ...
        ~isscalar(maxiter) || ~isfinite(maxiter) || maxiter < 0 || maxiter ~= floor(maxiter)
    error('Invalid fixed-point stopping parameters.');
end
history = zeros(maxiter+1,4);  % k, actual jump, update block norm, ratio
previous = NaN;
first = S;
converged = false;
reason = 'maximum_iterations';
for k = 0:maxiter
    [r,jump] = actual_jump_coordinates(S);
    ratio = NaN;
    if k > 0 && previous > 0 
        ratio = jump/previous; 
    end
    history(k+1,:) = [k,jump,0,ratio];
    if ~all(isfinite(r)) || ~isfinite(jump) 
        reason = 'nonfinite_jump'; 
        break; 
    end
    if jump <= tol
        converged = true; 
        reason = 'actual_jump_tolerance'; 
        break;
    end
    if k == maxiter 
        break; 
    end
    delta = solve_velocity_system(S,r/2);
    step = max([0;sqrt(sum(reshape(delta,S.M.dim,[]).^2,1))']);
    history(k+1,3) = step;
    if ~all(isfinite(delta)) 
        reason = 'nonfinite_update'; 
        break; 
    end
    previous = jump;
    S = set_spline_coefficients(S,S.c+delta);
    if k == 0 
        first = S; 
    end
end
S.maximum_actual_jump = jump;
S.fixed_point_converged = converged;
info = struct('converged',converged,'reason',reason,'iterations',k, ...
    'maximum_jump',jump,'tolerance',tol,'history',history(1:k+1,:));
end

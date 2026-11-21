function [x,v,a,b] = reference_evaluate(R,t,side)
%REFERENCE_EVALUATE Evaluate the piecewise-IVP reference in physical time.
% The two-sided position, velocity and acceleration agree at old knots.
% SIDE='left' requests the incoming jerk at a knot; default is 'right'.
if nargin < 3
    side = 'right';
end
tol = 2e-12*(1+R.T);
if t < -tol || t > R.T+tol
    error('Reference time outside its domain.');
end
t = max(0,min(R.T,t));
if strcmp(R.mode,'saved')
    D = R.dense;
    k = round(t/R.T*(numel(D.t)-1))+1;
    if abs(D.t(k)-t) > tol
        error('Saved replay uses configured grids; integrate for another grid.');
    end
    x = D.X{k};
    v = D.V{k};
    a = D.A{k};
    b = D.B{k};
    knot = round(t/R.H);
    if strcmp(side,'left') && knot >= 1 && abs(t-knot*R.H) < tol
        b = R.jerk_end{knot};
    end
    return;
end
ratio = t/R.H;
if abs(ratio-round(ratio)) < tol/R.H
    ratio = round(ratio);
end
j = min(R.n,max(1,floor(ratio)+1));
if strcmp(side,'left') && ratio == round(ratio) && ratio > 0
    j = round(ratio);
end
[x,v,a,b] = reference_state(R.M,deval(R.segments{j},t-R.t(j)));
end

function [z,dz,ddz] = euclidean_eval(P,t)
%EUCLIDEAN_EVAL Position and first two physical-time derivatives of a spline.
% Taking the interval on the right at a knot is harmless for this C2 spline.
j = find(P.t<=t,1,'last'); 
j = min(j,numel(P.t)-1);
if isempty(j) 
    error('Time precedes interpolation interval.'); 
end
h = P.t(j+1)-P.t(j); 
s = (t-P.t(j))/h;
p = P.X{j}; 
q = P.X{j+1}; 
v = P.V{j}; 
w = P.V{j+1};
z = (1-3*s^2+2*s^3)*p+(3*s^2-2*s^3)*q ...
    +h*(s-2*s^2+s^3)*v+h*(s^3-s^2)*w;
dz = (-6*s+6*s^2)*p/h+(6*s-6*s^2)*q/h ...
    +(1-4*s+3*s^2)*v+(3*s^2-2*s)*w;
ddz = (-6+12*s)*p/h^2+(6-12*s)*q/h^2 ...
    +(-4+6*s)*v/h+(6*s-2)*w/h;
end

function [x,v,a] = manifold_hermite(M,p,q,v0,v1,dt,s)
%MANIFOLD_HERMITE Left-centred exponential image of cubic Hermite data.
% s lies in [0,1]; v and a use physical time, not the normalized parameter.
% w=(d log_p)_q v1 is required for exact endpoint velocity interpolation.
% Factors dt^(-1) and dt^(-2) below convert s-derivatives into time derivatives.
u = M.log(p,q); 
w = M.dlog(p,q,v1);
z = (3*s^2-2*s^3)*u+dt*(s-2*s^2+s^3)*v0+dt*(s^3-s^2)*w;
dz = (6*s-6*s^2)*u/dt+(1-4*s+3*s^2)*v0+(3*s^2-2*s)*w;
ddz = (6-12*s)*u/dt^2+(-4+6*s)*v0/dt+(6*s-2)*w/dt;
[x,v,a] = M.jet(p,z,dz,ddz);
end

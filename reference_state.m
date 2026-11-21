function [x,v,a,b] = reference_state(M,y)
%REFERENCE_STATE Unpack one IVP state and remove numerical constraint drift.
% a=nabla_t v and b=nabla_t a are intrinsic vectors, not ordinary x''/x'''.
Y = reshape(y,numel(y)/4,4);
if strcmp(M.name,'sphere')
    x = Y(:,1)/norm(Y(:,1));
    v = Y(:,2)-x*(x'*Y(:,2));
    a = Y(:,3)-x*(x'*Y(:,3));
    b = Y(:,4)-x*(x'*Y(:,4));
else
    x = reshape(Y(:,1),M.shape);
    x = (x+x')/2;
    v = reshape(Y(:,2),M.shape);
    v = (v+v')/2;
    a = reshape(Y(:,3),M.shape);
    a = (a+a')/2;
    b = reshape(Y(:,4),M.shape);
    b = (b+b')/2;
end
end

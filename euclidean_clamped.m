function P = euclidean_clamped(t,X,v0,vn)
%EUCLIDEAN_CLAMPED Classical clamped cubic spline for vector/matrix cells.
% Solve for unknown interior velocities by equality of second derivatives.
% This is used for the initial-chart and ambient radial-projection baselines.
t = t(:); 
n = numel(t)-1; 
dt = diff(t); 
V = cell(n+1,1);
V{1} = v0; 
V{n+1} = vn;
if n > 1
    K = zeros(n-1); 
    rhs = zeros(n-1,numel(X{1}));
    for j = 2:n
        row = j-1;
        r = 3*(X{j+1}-X{j})/dt(j)^2+3*(X{j}-X{j-1})/dt(j-1)^2;
        K(row,row) = 2*(1/dt(j-1)+1/dt(j));
        if j == 2 
            r = r-v0/dt(j-1); 
        else
            K(row,row-1) = 1/dt(j-1); 
        end
        if j == n 
            r = r-vn/dt(j); 
        else
            K(row,row+1) = 1/dt(j); 
        end
        rhs(row,:) = r(:)';
    end
    W = K\rhs;
    for j = 2:n 
        V{j} = reshape(W(j-1,:),size(X{j})); 
    end
end
P = struct('t',t,'X',{X},'V',{V});
end

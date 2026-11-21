function E = discrete_cubic_energy(M,points,T)
%DISCRETE_CUBIC_ENERGY Scaled sum of centered logarithmic differences.
% points is a cell array sampled at K+1 UNIFORM physical times in [0,T].
% E = sum_i ||log_{p_i}p_{i-1}+log_{p_i}p_{i+1}||_{p_i}^2/(2*delta^3),
% delta=T/K; i runs over the interior energy samples only.
% The minus-time logarithm is added because it points backwards. The sum
% approximates delta^2 times covariant acceleration. Squaring and summing
% therefore requires delta^(-3), including the 1/2 from the cubic functional.
% Use the same K and T for every comparator. No acceleration jet is used.
K = numel(points)-1;
if K < 2 || T <= 0 
    error('Need positive duration and at least two subintervals.'); 
end
delta = T/K; 
values = zeros(K-1,1);
for i = 2:K
    a = M.log(points{i},points{i-1})+M.log(points{i},points{i+1});
    values(i-1) = M.inner(points{i},a,a);
end
E = sum(values)/(2*delta^3);
end

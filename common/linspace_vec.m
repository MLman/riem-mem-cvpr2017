function V = linspace_vec(a,b,n)
% LINSPACE_VEC generates n vectors between column vector a and b.
ndim = size(a,1);
V = zeros(ndim,n);
for i = 1:ndim
    V(i,:) = linspace(a(i),b(i),n);
end
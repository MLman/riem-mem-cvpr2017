function V = logmap_pt2array_spd(p,X)
%LOGMAP_PT2ARRAY_SPD returns logmap(p,Y) for SPD manifolds. This is the
%faster version of LOGMAP_VECS_SPD for one p not a set of base points P.
%
%    V = LOGMAP_PT2ARRAY_SPD(p,X)
%
%    p is a SPD matrix. (dim_p x dim_p, where dim_p = size(p,1))
%    X is a set of SPD matrices (dim_p x dim_p x N).
%    V is a set of symmetric matrices (dim_p x dim_p x N).
%
%   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD

%   $ Hyunwoo J. Kim $  $ 2016/10/25 21:28:52 (CDT) $
%   Warning is added.

[U, D] = eig(p);
g = U*sqrt(D);
invg = diag(1./sqrt(diag(D)))*U'; % 1.3 X faster

V = zeros(size(X));
%% For each data
for i = 1:size(X,3)
    if norm(p-X(:,:,i)) < 1e-18
        V(:,:,i) = zeros(size(p));
        continue
    end
    y = invg*X(:,:,i)*invg';
    [U, S] = eig(y);
    H = g*U;
    V(:,:,i) = H*diag(log(diag(S)))*H';
    V(:,:,i) = (V(:,:,i)+V(:,:,i)')/2; % make it symmetric.
end
if ~isfinite(V) | ~isreal(V)
    warning(sprintf('%s: Numerical Error. Nonreal tangent vectors.',mfilename));
end
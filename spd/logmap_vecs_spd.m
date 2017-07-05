function V = logmap_vecs_spd(X,Y)
%LOGMAP_VECS_SPD returns logmap(X,Y) for SPD manifolds.
%
%    V = LOGMAP_VECS_SPD(P,X)
%
%    X, Y is a set of SPD matrices (dimX x dimX x N, where dimX = size(X,1)).
%    V is a set of symmetric matrices.
%
%   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:42:13 $


V = zeros(size(Y));
if size(X,3) ==1 
    for i = 1:size(Y,3)
        yi = Y(:,:,i);
        V(:,:,i) = logmap_spd(X,yi);
    end
else
    for i = 1:size(X,3)
        xi = X(:,:,i);
        yi = Y(:,:,i);
        V(:,:,i) = logmap_spd(xi,yi);
    end
end


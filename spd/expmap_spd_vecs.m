function exp_p_x = expmap_spd_vecs(P,X)
%EXPMAP_SPD_VECS maps tangent vector X onto SPD manifold. X is a mxstack.
%
%    exp_p_x = EXPMAP_SPD(P,X)
%
%    P, exp_p_x is a SPD matrix.
%    X is a symmetric matrix.
%    
%   See also LOGMAP_SPD, KARCHER_MEAN_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:26:13 $

    [U,D] = eig(P);
    g = U*sqrt(D);
    invg = inv(g);
    exp_p_x = zeros(size(X));
    for i=1:size(X,3)
        if norm(X(:,:,i)) < 1e-18
            exp_p_x = P;
            return
        end
        Y = invg*X(:,:,i)*invg';
        [V, S] = eig(Y);
        gv = g*V;
        exp_p_x(:,:,i) = gv*diag(exp(diag(S)))*gv';
    end
end
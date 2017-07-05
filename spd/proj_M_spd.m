function p = proj_M_spd(X,varargin)
%PROJ_M_SPD projects a matrix onto SPD manifolds.
%
%
%    Example:
%        p = PROJ_M_SPD(X)
%
%   p is the point on SPD manifolds.
%   X is a n x n matrix.
%
%   See also MGLM_LOGEUC_SPD, MGLM_SPD

%   $ Hyunwoo J. Kim $  $ 2016/10/25 14:24:26 (CDT) $

    if nargin == 2
        c = varargin{1};
    else 
        c = eps;
    end
    
    % Make a matrix symmetric positive definite.
    if norm(X-X') > eps
        X = (X+X')/2;
    end
    [V, D ] = eig(X);
    D = diag(D);  
    if prod(D > c) % If all eig values are big enough, just return it.
        p = X;
        return
    end
    p = zeros(size(X));
    for i =1:length(D)
        if D(i) > 0
            p = p + D(i)*V(:,i)*V(:,i)';
        end
    end
    % Now X is spd
    % Make psd matrix
    if sum(D > 0+c) < length(D)
        a = 1e-16; 
        pnew = p;
        while ~isspd(pnew, c)
            pnew = p + a*eye(3);
            a = 2*a;
        end
        p = pnew;
    end
end
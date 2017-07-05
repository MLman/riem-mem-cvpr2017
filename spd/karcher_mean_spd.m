function xbar = karcher_mean_spd(X, W, niter)
%KARCHER_MEAN_SPD calculates the intrinsic mean with weight W on SPD manifolds.
%
%   xbar = KARCHER_MEAN_SPD(X, [], niter)
%   xbar = KARCHER_MEAN_SPD(X, W, niter)
%
%   W is weights.
%   X is a set of points on SPD manifolds.
%   xbar is the Karcher mean of X.
%   niter is the maximum iterations.
%
%   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD

%   $ Hyunwoo J. Kim $  $ 2016/10/25 20:57:04 (CDT) $
%   ISSPD is added.

% TODO: check the function value and find the step size to decrease it
% monotonically.

if isempty(W)
    xbar = mean(X,3);
    for iter = 1:niter
        phi = mean(logmap_pt2array_spd(xbar,X),3);
        % TODO: check xbar is SPD.
        xbar = expmap_spd(xbar, phi);
        assert(isspd(xbar)==1);
        if norm(phi) < 1e-18
            break
        end
    end
else
    [~, idx] = max(W);
    xbar = X(:,:,idx);
    W = W/norm(W,1);
    for iter = 1:niter
        tmp = logmap_pt2array_spd(xbar,X);
        wtmp = zeros(size(tmp));
        for i = 1:size(tmp,3)
            wtmp(:,:,i) = W(i)*tmp(:,:,i);
        end
        phi = sum(wtmp,3);
        xbar = expmap_spd(xbar, phi);
        assert(isspd(xbar)==1);
        if norm(phi) < 1e-18
            break
        end
    end
end

function Gr = uni_interpolation(A,B,r, varargin)
% UNI_INTERPOLATION performs univariate interpolation between A, and B 
% (0 <= r <= 1), with respect to a Riemannian metric. This interpolation is
% the same as the geometric mean w.r.t the Euclidean metric.
%
%   GM(r) = A*(inv(A)*B)^r = B*(inv(B)*A)^(1-r)
%   GM(0) == A and GM(1) == B
%
%   See Also: KARCHER_MEAN_SPD

%   $ Hyunwoo J. Kim $  $ 2016/10/25 11:53:29 (CDT) $
    assert(0 <= r && r <=1);
    if nargin >3 
        c = varargin{1};
    else
        c = 1e-15;
    end
    if r < 1
        Gr = A*(A\B)^r; % Make this SPD.
    else
        Gr = B*(B\A)^(1-r);
    end
    Gr = proj_M_spd(Gr,c);
%    Gr = B*(inv(B)*A)^(1-r);
end















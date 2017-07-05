function r = norm_TpM_spd(p,v)
%NORM_TPM_SPD calculates the norm of tangent vector v in TpM on SPD manifolds.
%
%    r  = NORM_TPM_SPD(p,v)
%
%   See also DIST_M_SPD, NORM_TPM_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:24:38 $
    r = sqrt(innerprod_TpM_spd(v,v,p));
end

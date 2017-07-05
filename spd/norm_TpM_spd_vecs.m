function r = norm_TpM_spd_vecs(p,vs)
%NORM_TPM_SPD_VECS evaluates norm of vs in TpM.
%
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2016/11/20 00:29:10 (CST) $
    r = sqrt(innerprod_TpM_spd_vecs(vs, vs, p));
end

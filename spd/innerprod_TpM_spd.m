function r = innerprod_TpM_spd(U,V,P)
%INNERPROD_TPM_SPD calculates the inner product of U and V in T_{P}M on SPD manifolds.
%
%    r  = INNERPROD_TPM_SPD(U,V,P)
%
%   See also DIST_M_SPD, NORM_TPM_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:23:38 $

    try
        invP = inv(P);
    catch
        invP = pinv(P);
        disp('pinv');
    end
    sqrtinvP= sqrtm(invP);
    r = trace(sqrtinvP*U*invP*V*sqrtinvP);
end
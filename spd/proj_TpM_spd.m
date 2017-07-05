function V = proj_TpM_spd(V)
%PROJ_TPM_SPD projects a set of tangent V vectors onto TpM. Symmetrization.
%
%   See also MGLM_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:59:20 $

    for i = 1:size(V,3)
        V(:,:,i) = (V(:,:,i)+V(:,:,i)')/2;
    end
end
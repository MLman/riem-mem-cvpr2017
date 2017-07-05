function W = innerprod_TpM_spd_wrapper(Y, V, P)
%INNERPROD_TPM_SPD_WRAPPER projects Y to the submanifold, EXP(P,V) by
% the innerproduct. This is useful for PGA type of visualization.
%
%
%   See Also: INNERPROD_TPM_SPD

%   $ Hyunwoo J. Kim $  $ 2016/10/25 20:40:18 (CDT) $

Ywr = logmap_vecs_spd(P, Y);
W = zeros(size(Y,3), size(V,3));
for i=1:size(Ywr,3)
    for j =1:size(V,3)
        W(i,j) = innerprod_TpM_spd(Ywr(:,:,i),V(:,:,j),P);
    end
end
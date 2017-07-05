function Vnew = embeddingRd_vecs_atI(V)
%EMBEDDINGRD_VECS_ATI
%
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2016/11/03 14:06:33 (CDT) $

    Vnew = symmx2vec(V);

    % Multiply sqrt(2) to off diagonal elements.
    dim = size(V,1);
    w = get_weights_for_symmx2vec(dim);    
    Vnew = repmat(w, 1, size(V,3)).*Vnew;
end
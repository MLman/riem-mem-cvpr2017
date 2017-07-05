function Vnew = embeddingRd_vecs(p,V)
%EMBEDDINGR6_VECS embeds a set of tangnet vectors V in TpM onto R6. This transports V from TpM to T_{I}M.
%
%   This will be replaced later with the group action for n x n spd matrix and a vectorization. 
%    
%   Vnew = EMBEDDINGRD_VECS(p,V)
%
%   See also EMBEDDINGR6_VECS, LOGMAP_SPD

%   $ Hyunwoo J. Kim $  $ 2016/11/03 13:41:12 (CDT) $
%   Weights are added.

    VwrI = grpaction_p2i(p, V);
    Vnew = symmx2vec(VwrI);

    % Multiply sqrt(2) to off diagonal elements.
    dim = size(p,1);
    w = get_weights_for_symmx2vec(dim);    
    Vnew = repmat(w, 1, size(V,3)).*Vnew;
end
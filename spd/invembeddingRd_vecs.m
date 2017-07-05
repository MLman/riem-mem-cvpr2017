function Vnew = invembeddingRd_vecs(p,Vecs)
%INVEMBEDDINGRD_VECS
%
%    
%
%
%   See Also: 

%   $ Hyunwoo J. Kim $  $ 2016/11/03 22:59:36 (CDT) $

%   $ Revision: 0.13 $  

    VatI = vec2symmx(Vecs);
    Vnew = grpaction_i2p(p, VatI);
    
    % Multiply sqrt(2) to off diagonal elements.
    dim = size(Vecs,1);
    nrows = (-1+sqrt(1+8*dim))/2;
    w = get_weights_for_vec2symmx(nrows);    
    Vnew = repmat(w, [1, 1,size(Vecs,2)]).*Vnew;
end
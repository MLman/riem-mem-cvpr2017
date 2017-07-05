function Vnew = invembeddingRd_vecs_atI(Vecs)
%INVEMBEDDINGRD_VECS_ATI is the inversion of EMBEDDINGRD_VECS_ATI.
%   
%   This will be replaced later with the group action for 
%   n x n spd matrix and a vectorization. 
%    
%   Vnew = INVEMBEDDINGRD_VECS_ATI(p,V)
%
%   See also EMBEDDINGRD_VECS_ATI

%   $ Revision: 0.13 $  
    VatI = vec2symmx(Vecs);
   
    % Multiply sqrt(2) to off diagonal elements.
    dim = size(Vecs,1);
    nrows = (-1+sqrt(1+8*dim))/2;
    w = get_weights_for_vec2symmx(nrows);    
    Vnew = repmat(w, 1, size(Vecs,3)).*VatI;
end
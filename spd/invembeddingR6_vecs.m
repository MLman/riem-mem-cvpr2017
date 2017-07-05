function Vnew = invembeddingR6_vecs(p, V)
%INVEMBEDDINGR6_VECS is the inverse of EMBEDDINGR6_VECS. This brings V in R6 to in TpM.
%This transports V from T_{I}M to TpM.
%
%   This will be replaced later with the group action for n x n spd matrix and a proper vectorization. 
%    
%   Vnew = EMBEDDINGR6_VECS(p,V)
%
%   p is a base point. 
%   Vnew is a set of tangent vectors in TpM.
%   V is a set of vectors in R6.
%
%   See also INVEMBEDDINGR6, EMBEDDINGR6_VECS, EMBEDDINGR6

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:32:53$ 

%    v = [Sxx Sxy Sxz Syy Syz Szz]';    

if size(V,2) == 1
    Vnew = invembeddingR6(p, V);
    return;
elseif length(size(V)) == 2
    nmx = size(V, 2);
    Vnew = zeros(3,3, nmx);
    for i=1:nmx
        Vnew(:,:,i) = invembeddingR6(p, V(:,i));
    end
else
    error('V is wrong input');
end
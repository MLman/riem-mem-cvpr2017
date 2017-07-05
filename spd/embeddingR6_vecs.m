function Vnew = embeddingR6_vecs(p,V)
%EMBEDDINGR6_VECS embeds a set of tangnet vectors V in TpM onto R6. This transports V from TpM to T_{I}M.
%
%   This will be replaced later with the group action for n x n spd matrix and a vectorization. 
%    
%   Vnew = EMBEDDINGR6_VECS(p,V)
%
%   See also INVEMBEDDINGR6_VECS, LOGMAP_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:16:53 

if length(size(V)) == 2
    Vnew = embeddingR6(p, V);
    return;
elseif length(size(V)) == 3
    nmx = size(V,3);
    Vnew = zeros(6,nmx);
    for i=1:nmx
        Vnew(:,i) = embeddingR6(p, V(:,:,i));
    end
else
    error('V is wrong input.');
end

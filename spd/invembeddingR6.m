function vnew = invembeddingR6(p, v)
%INVEMBEDDINGR6 is the inverse of EMBEDDINGR6. This brings v in R6 to in TpM.
%This transports v from T_{I}M to TpM.
%
%   This will be replaced later with the group action for n x n spd matrix and a proper vectorization. 
%    
%   vnew = EMBEDDINGR6(p,v)
%
%   p is a base point. 
%   vnew is a tangent vector in TpM.
%   v is a vector in R6.
%
%   See also INVEMBEDDINGR6_VECS, EMBEDDINGR6_VECS, EMBEDDINGR6

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:34:51$ 

%    v = [Sxx Sxy Sxz Syy Syz Szz]';    

%   step 1
    w = [1  sqrt(2) sqrt(2)   1 sqrt(2) 1]';
    S = vec2symmx(v./w);
    
%   step 2    
    sqrtp= sqrtm(p);
    vnew = sqrtp*S*sqrtp;
end
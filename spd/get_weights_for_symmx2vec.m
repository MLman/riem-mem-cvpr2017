function w = get_weights_for_symmx2vec(p)    
%GET_WEIGHTS_FOR_SYMMX2VEC
%
%
%   See Also: GET_WEIGHTS_FOR_VEC2SYMMX

%   $ Hyunwoo J. Kim $  $ 2016/11/03 13:47:14 (CDT) $

W = sqrt(2)*ones(p);
W = W - diag(diag(W))+eye(p); 
w = symmx2vec(W);
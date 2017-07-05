function W = get_weights_for_vec2symmx(p)    
%GET_WEIGHTS_FOR_VEC2SYMMX
%
%
%   See Also: GET_WEIGHTS_FOR_SYMMX2VEC

%   $ Hyunwoo J. Kim $  $ 2016/11/03 13:48:52 (CDT) $

    W = 1./sqrt(2)*ones(p);
    W = W - diag(diag(W))+eye(p); 
end
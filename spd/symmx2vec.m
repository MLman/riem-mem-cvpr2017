function v = symmx2vec(mx)
%SYMMX2VEC converts matrices mx to vectors v. 
%
%   v = symmx2vec(mx)
%
%   v is a set of n(n+1)/2 dimensional vectors to n by n matrices.
%   mx is a set of n x n matrices.
%
%   See also INVEMBEDDINGR6, VEC2SYMMX

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:09:53 $

    [ nrow, ncol, ndata ] = size(mx);
    v = zeros(nrow*(nrow+1)/2,ndata);
    k =1;
    for i=1:ncol
        for j=i:ncol
            v(k,:) = squeeze(mx(i,j,:))';
            k = k + 1;
        end

    end
end
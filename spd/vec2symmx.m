function mx = vec2symmx(v)
%VEC2SYMMX converts vectors v to matrices mx.
%
%   mx = vec2symmx(v)
%
%   v is a set of n(n+1)/2 dimensional vectors to n by n matrices.
%   mx is a set of n x n matrices.
%
%   See also INVEMBEDDINGR6, SYMMX2VEC

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:09:53 $

    [dimv ndata] = size(v);
    n = (-1 + sqrt(1+8*dimv))/2;
    mx = zeros(n,n,ndata);
    k = 1;
    for i=1:n
        for j=i:n
            mx(i,j,:) = v(k,:);
            if i ~=j
                mx(j,i,:) = v(k,:);
            end
            k = k + 1;
        end
    end
end
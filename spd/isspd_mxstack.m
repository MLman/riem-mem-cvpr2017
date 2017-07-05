function t = isspd_mxstack(Y,varargin)
t = 0 ;
if nargin < 2
    for i=1:size(Y,3)
        t =t + (~isspd(Y(:,:,i)));
    end
else
    c = varargin{1};
    for i=1:size(Y,3)
        t =t + (~isspd(Y(:,:,i),c));
    end
end
t = ~t;
function Yhat = prediction_spd(p,V,X)
%PREDICTION_SPD predicts phat based on estimate p, V and covariate X.
%
%   p is a base point (SPD maxtrix). 
%   V is a set of tangent vectors (3 x 3 x dimX symmetric matrix).
%   X is a set of covariates, dimX x N column vectors.
%   p_hat is the prediction.
%
%   See also MGLM_SPD, FEVAL_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 00:13:20 $
    
[ndimX ndata] = size(X);
Yhat = zeros(size(p,1),size(p,2),ndata);

for i = 1:ndata
    Vi = zeros(size(p));
    for j = 1:ndimX
        Vi = Vi+V(:,:,j)*X(j,i);
    end
    Yhat(:,:,i) = expmap_spd(p,Vi);
end

function E = feval_spd(p,V,X,Y)
%FEVAL_SPD evaluates the objective function value (the sum of squared geodesic errors) of MGLM on SPD. 
%
%    E = FEVAL_SPD(p,V, X, Y)
%
%   !! make sure that X is centered if p, V are calculated by centered X !!
%
%   X is a set of column vectors. (dimX X N, where dimX = size(X,1))
%   p is the base point.
%   V is a set of symmetric matrices (dim_p x dim_p x dimX).
%   Y is a set of SPD matrices (dim_p x dim_p x N, where dim_p = size(p,1)).
%   E is the sum of squared geodesic errors.
%
%   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD, GSQERR_SPD,
%   PREDICTION_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 15:42:13 $

P_hat = prediction_spd(p,V,X);
E = gsqerr_spd(Y, P_hat);

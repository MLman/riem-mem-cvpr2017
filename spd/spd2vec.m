function v = spd2vec(X)
%SPD2COEF converts matrix X to coefficient vector v.
%
%   v = [Dxx, Dyx, Dyy, Dzx, Dzy, Dzz]' 
%
%   See Also: DTICOEF2MX
%
%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:56:53 $

Dxx = X(1,1);
Dyx = X(2,1);
Dyy = X(2,2);
Dzx = X(3,1);
Dzy = X(3,2);
Dzz = X(3,3);
v = [Dxx, Dyx, Dyy, Dzx, Dzy, Dzz]';

function r2  = r2stat_spd(Y_bar, Y, Y_hat)
%R2STAT_SPD calculates R2 statistic.
%
%    r2  = R2STAT_SPD(Y_bar, Y, Y_hat)
%
%    Example:
%        Y, Y_hat; % Given Y, Y_hat
%        Y_bar = karcher_mean_spd(Y,[],1000);
%        r2  = r2stat_spd(Y_bar, Y, Y_hat);
%
%   See also GSQERR_SPD, FEVAL_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:03:38 $

gvar = gsqerr_spd(repmat(Y_bar,[1,1,size(Y,3)]), Y);
uvar = gsqerr_spd(Y, Y_hat);
r2 = 1-uvar/gvar;

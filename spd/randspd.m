function P = randspd(n, varargin)
%RANDSPD generates n by b random symmatrix positive matrix P.
%
%   P = randspd(n)
%   P = randspd(n,c)
%   P = randspd(n,c,udist)
%
%   c is parameter for variance. Bigger c has bigger variance.
%   udist is the upper bound of distance from I to P w.r.t. GL-invariant
%   measure.
%
%   See also SYNTH_DTI_DATA

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:03:38 $
if nargin >= 2
    c = varargin{1};
else
    c = 3;
end

P = c*(rand(n)-0.5);
P = P*P';

if nargin == 3
    udist = varargin{2};
    while dist_M_spd(P,eye(n)) > udist
        P = c*(rand(n)-0.5);
        P = P*P';
    end
end

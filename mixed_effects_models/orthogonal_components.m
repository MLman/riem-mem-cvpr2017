function VV = orthogonal_components(p, u, V)
%ORTHOGONAL_COMPONENTS get orthogonal components
% Assuming U is a normal tangnet vector at P. This returns a component of V 
% which is perpendicular to u. 
% 
%
%   See Also: 

%   $ Hyunwoo J. Kim $  $ 2016/11/02 23:53:21 (CDT) $

a = innerprod_TpM_spd(u, V, p);
VV = V- a*u;
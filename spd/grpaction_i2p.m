function Vp = grpaction_i2p(p, Vs)
%GRPACTION_I2P performs group actions for each V of Vs to translate from
%T_{I}M to T_{p}M. p^1/2 V p^T/2 = p^1/2 V p^1/2 since p^1/2 is symmetric.
%
%
%   See Also: GRPACTION_P2I,QUADFUNC

%   $ Hyunwoo J. Kim $  $ 2016/04/20 15:06:27 (CDT) $
    [U,S] = eig(p);
    pp = U*diag(sqrt(diag(S)))*U';
    Vp = quadfunc(pp, Vs);
end
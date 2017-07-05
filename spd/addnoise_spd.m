function Anew = addnoise_spd(A, maxerr)
V = randsym(size(A,1));
if norm_TpM_spd(A,V) > maxerr
    V = V/norm_TpM_spd(A,V)*maxerr;
end
Anew = expmap_spd(A, V);

function M = randsym(n)
M = randn(n);
M = (M+M')/2;
function safev = zero2one(v)
% ZERO2ONE is a safegauard function for division by zero. 
if v == 0
    safev = 1;
else
    safev = v;
end

function T = ismonoinc(v)
% Returns whether T is monotonically increasing or not. 
% This rules out constant functions.
% 
if length(v) == 1 % Constant function
    T = false;
    return
end

incs = v(2:end)-v(1:end-1);
if norm(incs) == 0 || sum(incs < 0) > 0
 % constant function OR decreasing 
    T = false;
    return 
end

T = true;
return  

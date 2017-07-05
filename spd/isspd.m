function T = isspd(mx,varargin)
%ISSPD check mx is a symmetric positive definite matrix.
%    This check whether the smallest eigen value is bigger than c.
%    Default c is epsilon.
%
%    Example:
%        T = isspd(mx)
%        T = isspd(mx,C)
%
%   See also MGLM_LOGEUC_SPD, MGLM_SPD, PROJ_M_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 16:40:17$ 

    if nargin ==2
        c = varargin{1};
    else 
        c = eps;
    end
    % Check matrices are symmetric positive definite.
    T = zeros(size(mx,3),1);
    for i=1:size(mx,3)
        T(i) = (sum(eig(mx(:,:,i)) <= 0+c ) ==0) && issym(mx(:,:,i));
    end
end

function [T S] = issym(mx)
    tol = 0.00001;
    S = zeros(size(mx,3),1);
    for i = 1:size(mx,3)
        S(i) = (sum(sum(abs(mx(:,:,i)-mx(:,:,i)'))) < tol);
    end
    T = (sum(S) == size(mx,3));
end
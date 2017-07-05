function P = myrandperm(N, K)
%MYRANDPERM returns N+1 different permutations including 1:K at the top of matrix P.
%
%   P in R^(N+1)xK
%
%   See Also: RANDPERM

%   $ Hyunwoo J. Kim $  $ 2016/08/13 16:21:05 (CDT) $
%   $ Hyunwoo J. Kim $  $ 2015/08/01 11:11:53 (CDT) $

%keyboard;
if(N>factorial(K))
    fprintf('Number of permutations isnt less than (num.subjects)! \nPlease add more subjects. Bye :-)\n');
    return
end
    P = zeros(N+1, K);
    % To be safe when randperm gives redundant permutations.
    num_perm = round(1.1*N); 
    P(1,:) = 1:K; %first perm is unpermuted
    for i=2:(num_perm+1)
        P(i,:)=randperm(K); %Guarantee of distinct permutations??
    end
    
    % Guarantee of distinct permutations, yes.
    % Due to unique, P looks ordered.
    P = unique(P,'rows'); 
    P = P(1:(N+1),:);
    
    % What is the below line of code doing?
    % Check the first row is not permuted.
    assert(sumabs(P(1,:) ~= 1:K) == 0);
end
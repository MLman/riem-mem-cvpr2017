function U = main_change_direction(X,Y,Z,option)
%MAIN_CHANGE_DIRECTION
%
%
%   See Also: MMEM

%   $ Hyunwoo J. Kim $  $ 2016/11/02 21:55:29 (CDT) $

centerX = optread(option, 'centerX', true);

nsubjects = size(Z,2);
nsamples = size(X,1);

if centerX
    Xbar = mean(X,1);
    Xci = zeros(nsubjects,1);
    for i=1:nsubjects
        Xci(i,:) = mean(X(Z(:,i)==1,:),1);
    end
    Xc = zeros(size(X));
    for i=1:nsamples
        Xc(i,:) = X(i,:) - Xci(Z(i,:)==1,:);
    end
    Xc = Xc';
else
    Xc = X';
end

dimY = size(Y,1);
%Ybar = karcher_mean_spd(Y,[],100);

Yibar = zeros(dimY, dimY, nsubjects);
for isub = 1: nsubjects
    Yis = Y(:,:,Z(:,isub) ==  1);
    Yibar(:, :, isub) = km_spd_wrapper(Yis, [], 50);
end

% Tanget vectors at pi(r)
Ywr = zeros(size(Y));
for i =1:nsamples
    Ywr(:,:,i) = logmap_spd(Yibar(:,:,Z(i,:)==1),Y(:,:,i));
end
I = eye(dimY);
% Group actions from pi(r) to p and then group action from pi(r) to I 
YwrAtP = zeros(size(Ywr));
for i =1:nsamples
    YwrAtP(:,:,i) = paralleltranslateAtoB_spd(Yibar(:,:,Z(i,:)==1), I, Ywr(:,:,i));
end
% Group action is NOT needed.
Yv = symmx2vec(YwrAtP);

% Perform MGLM with y, and Xij (or Xij(r)) 
L = Yv/Xc; % Peform MGLM. Centered.
U = vec2symmx(L);

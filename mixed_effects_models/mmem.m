function [Yhat, Bhat, P] = mmem(X,Y,Z, option)
%MMEM performs the Manifold-valued Mixed Effects Model.
%
%   Y_ij = EXP(pi, Gamma_{p->pi}(B)X_ij), (default)
%   Y_ij = EXP(pi, Gamma_{p->pi}(B) (X_ij- bar{X}i)) 
%   
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2016/10/25 22:49:13 (CDT) $
centerX = optread(option, 'centerX', true);
r = optread(option, 'r', 1);

nsubjects = size(Z,2);
nsamples = size(X,1);

if centerX
    Xbar = mean(X,1);
    Xci = zeros(nsubjects,1);
    for i=1:nsubjects
        Xci(i,:) = mean(X(Z(:,i)==1,:),1);
    end
    Xoff = zeros(size(X));
    for i = 1:nsamples
        Xoff(i,:) = (r*Xci(Z(i,:)==1,:)+(1-r)*Xbar);
    end
end



if centerX
    Xc = zeros(size(X));
    for i=1:nsamples
        Xc(i,:) = X(i,:)-Xoff(i,:);
    end
    Xc = Xc';
else
    Xc = X'; % Xc is the centered column vectors.
end

dimX = size(X,2);
dimY = size(Y,1);
Yibar = zeros(dimY, dimY, nsubjects);

for isub = 1: nsubjects
    Yis = Y(:,:,Z(:,isub) ==  1);
    Yibar(:, :, isub) = km_spd_wrapper(Yis, [], 50);
end

Ybar = karcher_mean_spd(Y,[],100);
P = Ybar;
% Calculate pi(r)

Pis = zeros(dimY,dimY, nsubjects);
for i=1:nsubjects
    Pis(:,:,i) = uni_interpolation(Ybar, Yibar(:,:,i),r);
end

% Tanget vectors at pi(r)
Ywr = zeros(size(Y));
for i =1:nsamples
    Ywr(:,:,i) = logmap_spd(Pis(:,:,Z(i,:)==1),Y(:,:,i));
end

% Group actions from pi(r) to p and then group action from pi(r) to I 
YwrAtP = zeros(size(Ywr));
for i =1:nsamples
    YwrAtP(:,:,i) = paralleltranslateAtoB_spd(Pis(:,:,Z(i,:)==1), P, Ywr(:,:,i));
end

% Group actions from P to I (?) if needed.
Yv = embeddingRd_vecs(P, YwrAtP);

% Perform MGLM with y, and Xij (or Xij(r)) 
L = Yv/Xc;% Peform MGLM. Centered.

logYv_hat = L*Xc; % Prediction.
% Transport V back to TpM
Bhat = invembeddingRd_vecs(P,L);

% Transport Y^{\wr} back to TpM
V_hat = invembeddingRd_vecs(P,logYv_hat);

% Transport Y^{\wr} to TpiM
Ywrhat = zeros(size(Ywr));
for i =1:nsamples
    Ywrhat(:,:,i) = paralleltranslateAtoB_spd(P, Pis(:,:,Z(i,:)==1), V_hat(:,:,i));
end

Yhat = zeros(size(Y));
for i =1:nsamples
    Yhat(:,:,i) = expmap_spd(Pis(:,:,Z(i,:)==1), Ywrhat(:,:,i));
end

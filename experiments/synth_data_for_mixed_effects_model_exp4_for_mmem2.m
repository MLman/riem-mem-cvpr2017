function [X, Y, Z, P, uV, B, Yintact, acc] = synth_data_for_mixed_effects_model_exp4_for_mmem2(option)
% Synthetic data with noise.
% For visualization, subjects are regularly spaced.
centerX = optread(option, 'centerX', true);
a1 = optread(option, 'a1', .5); % global slope 
a2 = optread(option, 'a2', -0.35); % local slope
a3 = optread(option, 'a3', 0.1); % std of local slopes
c1 = optread(option, 'c1', 2); % variance between subjects in X.
c2 = optread(option, 'c2', 0.8); % variance within a subject in X
c3 = optread(option, 'c3', 0.2); % noise

nsubjects = optread(option, 'nsubjects', 10);
nvisits = optread(option, 'nvisits', 4);
dimY = optread(option, 'dimY', 3);
npivots = optread(option, 'npivots', 2);
nsamples = nsubjects*nvisits;


% Longitudinal dependency
Z = zeros(nsamples, nsubjects); 
isub = 1;
nnvisit = 0; 
for i = 1:nsamples
    nnvisit = nnvisit + 1;
    Z(i,isub) = 1;
    if nnvisit >= nvisits
        nnvisit = 0; 
        isub = isub + 1;
    end
end

% Covariate generation
%Xci = c1*randn( nsubjects, 1);
Xci = c1*linspace(0,1,nsamples)';
X = zeros(nsamples,1);
for i=1:nsamples
    X(i,:) = Xci(Z(i,:)==1,:)+c2*randn;
end

% Pivot points
Yp = zeros(dimY, dimY, npivots+1);
Yp(:,:,1) = randspd(dimY,2,3);
for i=2:npivots
    Yp(:,:,i) = randspd(dimY,2,10);
end
P = Yp(:,:,1);

% Tangent vectors, geodesic bases for visualization.
V = zeros(3,3,npivots-1);
for j =1:npivots-1
    V(:,:,j) = logmap_spd(Yp(:,:,1),Yp(:,:,j+1));
end

B = V;

% Normalization of Bs, the main direction of change.
for i = 1:size(B,3)
    B(:,:,i) = B(:,:,i)/norm_TpM_spd(P, B(:,:,i));
end

u = randn(3,1); % Needed for many samples

Pi = zeros(dimY,dimY,nsubjects);
for i =1:nsubjects
    Pi(:,:,i) = expmap_spd(P, a1*B(:,:,1)*Xci(i,1));
end

r = 1; % Averaging portion, model assumption, 
       % r = 1; Each subject has its own intercept.
       % r = 0; All subjects have the same intercept.
       
%% Generate Ground Truth Data
B_news  = zeros(dimY, dimY, nsubjects);
acc = zeros(nsubjects,1);
for isub=1:nsubjects
    acc(isub) = (1+a3*randn(1));
    B_news(:,:,isub) = paralleltranslateAtoB_spd(P, Pi(:,:,isub), acc(isub)*a2*B(:,:,1));
end

Y = zeros(dimY,dimY,nsamples);
for i = 1:nsamples
    isub = find(Z(i,:));
    if centerX
        Y(:,:,i) = expmap_spd(Pi(:,:,isub), B_news(:,:,isub)...
            *(X(i,1)-Xci(Z(i,:)==1,1)));
    else
        Y(:,:,i) = expmap_spd(Pi(:,:,isub), B_news(:,:,isub)*X(i,1));
    end
end

%% Add noise
% TODO : Better sampling method is needed.
Yintact = Y; 
for i=1:nsamples
    Y(:,:,i) = addnoise_spd(Y(:,:,i),c3); % spdlognorm can be used. 
end

%% Sanity check
assert(isspd_mxstack(Y)==1, 'One of Y is not SPD');
uV = V;
for i = 1:size(V,3)
    uV(:,:,i) = uV(:,:,i) / norm_TpM_spd(P, uV(:,:,i));
end


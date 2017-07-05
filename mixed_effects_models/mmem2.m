function [Yhat, model] = mmem2(X, Y, Z, option)
%MMEM2 Riemannian Mixed Effects Models
%   
%
%   See Also:

%   $ Hyunwoo J. Kim $  $ 2016/11/02 22:53:21 (CDT) $
Z = logical(Z);
centerX = optread(option, 'centerX', true);
if centerX
    Xbar = mean(X,1);
    X = X-repmat(Xbar,size(X,1),1);
end

KMniter = optread(option, 'KMniter', 100);
niter = optread(option, 'niter', 100);
orthogonal_Ui = optread(option, 'orthogonal_Ui', true);

nsubjects = size(Z,2);
nsamples = size(X,1);
dimY = size(Y,1);

% STEP 1: Calculate the Population Frechet mean of Y. bar{y}
Ybar = karcher_mean_spd(Y, [], KMniter);
B=Ybar;

% STEP 2: Calculate the Frechet mean for each subject. bar{y}_i 
Yibar = zeros(dimY, dimY, nsubjects);
for isub = 1: nsubjects
    Yis = Y(:,:,Z(:,isub) ==  1);
    Yibar(:, :, isub) = km_spd_wrapper(Yis, [], 50);
end

% STEP 3: Main longitudinal change direction Eta by algorithm 3
% Eta in T_{I}M
Eta = main_change_direction(X,Y,Z,option);

% STEP 4: Calculate subject-specfic base points (random effects)
% B_i = EXP(B, U_i^{*})
% EtaAtB = Gamma_{I -> B}(eta)
% U_i  = argmin_{U_i \perp etaAtB } d( \bar{y}_i, EXP(B, U_i))^2

I = eye(dimY);
EtaAtB = paralleltranslateAtoB_spd(I, B, Eta);
EtaAtB = EtaAtB/norm_TpM_spd(B, EtaAtB); % Use only direction (unit vector).
Yibarwr = logmap_vecs_spd(B, Yibar);
Ui = zeros(size(Yibarwr));
for i =1:nsubjects
    if ~orthogonal_Ui
        Ui(:,:,i) = Yibarwr(:,:,i);
    else
        Ui(:,:,i) = orthogonal_components(B, EtaAtB, Yibarwr(:,:,i));
    end
end
Bi = zeros(size(Ui));
for i =1:nsubjects
    Bi(:,:,i) = expmap_spd(B, Ui(:,:,i));
end

%% Visualize here

% alppha = ones(nsubjects,1);
% tau = zeros(nsubjects, 1);
alppha = randn(nsubjects,1)*0.0001+1;
tau = randn(nsubjects,1)*0.0001;
% STEP 5: Tangent vectors to T_{I}M 
Ywr = zeros(size(Y));
for i =1:nsamples
    p = Bi(:,:,Z(i,:)==1);
    Ywr(:,:,i) = paralleltranslateAtoB_spd(p,I,logmap_spd(p, Y(:,:,i)));
end

% TODO: off diagonal vectors should be \sqrt{2}, since ..twices..
Yv = embeddingRd_vecs_atI(Ywr)';
dimYv = dimY*(dimY+1)/2;
Etav = embeddingRd_vecs_atI(Eta);

% For debugging, 
c = 1;
V = c*Etav;
t0 = 0; % It is undefined.


for iter=1:niter
    % STEP 8: Calculate teh subject-specific acceleration alpha_i and time
    % shift tau_i by generalized least square estimation with the priors
    % No priors for now.
%     le_error0=0;
%     for isam =1:nsamples
%         isub=find(Z(isam, :)==1);
%         le_error0 = le_error0+sum((Yv(isam,:)'- V*(alppha(isub)*(X(isam)-tau(isub)-t0)+t0)).^2);
%     end
    for isub = 1:nsubjects
        Yvi = Yv(Z(:,isub)==1,:);
        YYi = Yvi - repmat(V'*t0, size(Yvi,1),1);
        Wi = kron(X(Z(:,isub)), V')-t0*kron(ones(size(Yvi,1),1),V');
        Vi = kron(true(size(Yvi,1),1), V');
        % Should be generalized least square estimation
        A  =zeros(2);
        A(1,1) = sum(sum(Wi.*Wi));
        A(1,2) = -sum(sum(Wi.*Vi));
        A(2,1) = -A(1,2);
        A(2,2) = - sum(sum(Vi.*Vi));
        BB = [sum(sum(YYi.*Wi)); sum(sum(YYi.*Vi))];
        bb = A\BB;   
        alppha(isub) = bb(1);
        tau(isub) = bb(2)/bb(1);
    end
    %alppha
    %tau
    %fprintf('Iter %d\n',iter);
    
%    STEP 7: Calculate the common speed of change V=a\eta
%    and common time intercept t_0= b/(1-a).
    % Error check
%     le_error=0;
%     for isam =1:nsamples
%         isub=find(Z(isam, :)==1);
%         le_error = le_error+ sum((Yv(isam,:)'- V*(alppha(isub)*(X(isam)-tau(isub)-t0)+t0)).^2);
%     end
    Pij = zeros(dimYv,nsamples);
    Qi = zeros(dimYv, nsubjects);
    for isub=1:nsubjects
        %Qi(:,isub) = Etav *(1-alppha(isub));
        Qi(:,isub) = V *(1-alppha(isub));
    end
    
    for isam=1:nsamples
        isub = find(Z(isam, :)==1);
        %Pij(:,isam) = Etav*(alppha(isub)*(X(isam) - tau(isub)));
        Pij(:,isam) = V*(alppha(isub)*(X(isam) - tau(isub)));
    end
    
% |sum_{ij}q_i^Tq_i     & -sum_{ij}p_{ij}^Tq_i\\ | |b| = |sum_{ij} q_i^Ty^\wr_{ij}  |
% |-sum_{ij}p_{ij}^Tq_i & s um_{ij}p_{ij}^Tp_{ij}| |c|   |sum_{ij}p_{ij}^Ty^\wr_{ij}|

    A = zeros(2);
    BB = zeros(2,1); %The right hand side of the system of equations.
    qTq = zeros(nsubjects,1);
    for isub=1:nsubjects
       qTq(isub) = Qi(:,isub)'*Qi(:,isub);
    end
    nvisits = sum(Z,1)';
    A(1,1) = sum(qTq.*nvisits);
    pTq = 0;
    for isam = 1:nsamples
        pTq = pTq + Pij(:,isam)'*Qi(:,Z(isam,:)==1);
    end
    A(1,2) = pTq;
    A(2,1) = pTq;
    pTp = 0;
    for isam = 1:nsamples
        pTp = pTp + Pij(:,isam)'*Pij(:,isam);
    end
    A(2,2) = pTp;
    
    for isam =1:nsamples
        BB(1) = BB(1) + Qi(:,Z(isam,:)==1)'*Yv(isam,:)';
    end
    
    for isam =1:nsamples
        BB(2) = BB(2) + Pij(:,isam)'*Yv(isam,:)';
    end
    xx =  A\BB;
    b = xx(1);
    c = xx(2); 
    %fprintf('b=%f,c=%f, t0=%f\n',b, c,t0);

    V = c*Etav;
    if norm(c-1) < 1e-15
        %warning('Numerical Error');
        t0 = 0; % It is undefined.
    else
        t0 = b/c;
    end
%     le_error2=0;
%     for isam =1:nsamples
%         isub=find(Z(isam, :)==1);
%         le_error2 = le_error2+ sum((Yv(isam,:)'- V*(alppha(isub)*(X(isam)-tau(isub)-t0)+t0)).^2);
%     end
    %fprintf('le_error0=%f, le_error=%f, lerror2=%f\n',le_error0,le_error, le_error2);
end

model=[];
model.alppha = alppha;
model.tau = tau;
model.V = paralleltranslateAtoB_spd(I, B, invembeddingRd_vecs_atI(V));
model.B = B;
model.t0 = t0;
model.Ui = Ui;
model.centerX = centerX;
model.EtaatB = paralleltranslateAtoB_spd(I, B, Eta);
model.orthogonal_Ui = orthogonal_Ui;
Yhat = predict_mmem2(model, X, Y, Z);
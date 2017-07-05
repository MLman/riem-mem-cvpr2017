function X = spdlognorm(M, sig, varargin)
%SPDLOGNORM performs sampling on the tangent space at M. 
%   sig is the scalar value to control the deviation of samples.
%
%
%   See Also: spdlognormpdf

%   $ Hyunwoo J. Kim $  $ 2016/08/19 11:07:11 (CDT) $

NFAILS = 1000;
if nargin == 2
    nsamples = 1;
end
if nargin >= 3
    nsamples = varargin{1};
end
if nargin == 4
    udist = varargin{2}; %For truncated normal distribution.
end
dimM = (size(M,1)*(size(M,1)+1))/2;
if nargin <= 3
    v = mvnrnd(zeros(dimM,1),eye(dimM)*sig,nsamples)';    
    % Be careful with the order of elementers in v 
    % when covariance is not isotropic.
    V = vec2symmx(v); 
end
if nargin == 4
    V = zeros([size(M),nsamples]);
    for i =1:nsamples
        lVlp = inf;
        nfails = NFAILS;
        while lVlp > udist && nfails > 0
            v = mvnrnd(zeros(dimM,1),eye(dimM)*sig)';    
            Vi = vec2symmx(v);
            lVlp = norm_TpM_spd(M,Vi);
            nfails = nfails - 1;
        end
        assert(nfails > 0,'Too many rejections.');
        V(:,:,i) = Vi;
    end
end
X = expmap_pt2array_spd(M,V);
end
function w_new = paralleltranslateAtoB_spd(a, b, w)
%PARALLELTRANSLATEATOB_SPD transports a set of tangent vectors w from TaM to
%TbM.
%
%   w_new = PARALLELTRANSLATEATOB_SPD(a, b, w)
%
%   a, b are points on SPD matrices. 
%   w is a set of tangent vectors.
%   w_new is a set of transported tangent vectors.
%
%   See also MGLM_SPD

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 17:06:21 $

    if size(a,3) < size(b,3)
        a = repmat(a, [1 1 size(b,3)]);
    elseif size(a,3) > size(b,3)
        b = repmat(b, [1 1 size(a,3)]);
    end

    if size(b,3) ~= size(w,3)
        % a, b are fixed
        % This changes only w.
        fixab = 1;
        P1 = a;
        P2 = b;
    else
        fixab = 0;
    end

    w_new = zeros(size(w));
    if fixab
        if norm(a-b) < 1e-20
            w_new = w;
            return
        end
        rtp = sqrtm(a);
        invrtp = inv(rtp);
        v = logmap_spd(a,b);
        r = expm(invrtp*v/2*invrtp);
        L = rtp*r*invrtp;
        R = invrtp*r*rtp;
        for i= 1:size(w,3)
            w_new(:,:,i) = L*w(:,:,i)*R;
        end
    else
        for i = 1:size(w,3)
            P1 = a(:,:,i);
            P2 = b(:,:,i);
            if norm(P1-P2) < 1e-20
                w_new(:,:,i) = w(:,:,i);
                continue
            end
            w_new(:,:,i) = parallel(P1,P2,w(:,:,i));
        end
    end

    %% symmetrization.
    for i = 1:size(w,3)
        w_new(:,:,i) = (w_new(:,:,i)+w_new(:,:,i)')/2;
    end
end

function w_new = parallel(p,q,w)
    rtp = sqrtm(p);
    invrtp = inv(rtp);
    v = logmap_spd(p,q);
    r = expm(invrtp*v/2*invrtp);
    w_new = rtp*r*invrtp*w*invrtp*r*rtp;
end

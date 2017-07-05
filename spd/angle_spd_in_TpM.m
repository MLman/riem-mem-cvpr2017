function q = angle_spd_in_TpM(U,V,P)
lenU = norm_TpM_spd(P, U);
lenV = norm_TpM_spd(P, V);
UdotV = innerprod_TpM_spd(U, V, P);
% For numerical errors.
q = acos(max(min(UdotV/(lenU*lenV),1),-1));
function synth_plot_wrapper(X,Z,Y,uV,Pstar,P,Bhat,Yhat,P0,Bhat0,Yhat0,fname, mytitle)
W=innerprod_TpM_spd_wrapper(Y, uV, Pstar);
What=innerprod_TpM_spd_wrapper(Yhat,uV, Pstar);
What0=innerprod_TpM_spd_wrapper(Yhat0,uV, Pstar);
bhat=innerprod_TpM_spd(paralleltranslateAtoB_spd(P, Pstar,Bhat), uV, Pstar);
bhat0=innerprod_TpM_spd(paralleltranslateAtoB_spd(P0, Pstar,Bhat0), uV, Pstar);
synth_plot(X, W, Z, What, What0, bhat, fname,mytitle);
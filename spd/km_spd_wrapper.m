function xbar = km_spd_wrapper(X, W, niter)
    c = 1e-15;
    xbar = karcher_mean_spd(X, W, niter);
    xbar = proj_M_spd(xbar,c);
end
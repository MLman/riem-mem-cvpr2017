%function synth_exp2-2
clear
option=[];
option.c1=20;
option.c2=0.4;
option.c3=0.1;
option.a2 = -0.2;
[X, Y, Z, Pstar, uV, B] = synth_data_for_mixed_effects_model_exp4_for_mmem2(option);
W = innerprod_TpM_spd_wrapper(Y, uV, Pstar);


nsamples = size(X,1);
subjIDs = cell(nsamples,1);
for i =1:nsamples
    subjIDs{i} = num2str(find(Z(i,:) ==1));
end

figure
plot(X, W,'o');
hold on 
text(X, W, subjIDs,'VerticalAlignment','bottom','HorizontalAlignment','right')
title('Synthetic data.')
hold off

fname=sprintf('./synthdata/%s_%s.mat',mfilename,datestr(now,'yyyymmdd_HHMMSS'));
save(fname,'-v7.3');
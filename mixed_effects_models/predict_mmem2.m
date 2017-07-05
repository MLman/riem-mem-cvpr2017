function Yhat = predict_mmem2(model, X, Y, Z)
Ui = model.Ui;
B = model.B;
Bi = zeros(size(Ui));
t0 = model.t0;
alppha = model.alppha;
tau = model.tau;
[nsamples, nsubjects]=size(Z);
for i =1:nsubjects
    Bi(:,:,i) = expmap_spd(B, Ui(:,:,i));
end
V = model.V;
Vi = zeros(size(Bi));

for i = 1:nsubjects
    Vi(:,:,i) = paralleltranslateAtoB_spd(B, Bi(:,:,i), V);
end

T = Z*alppha .*(X-Z*tau-t0)+t0;
Yhat = zeros(size(Y));
for i =1:nsamples
    Yhat(:,:,i) = expmap_spd(Bi(:,:,Z(i,:)==1), Vi(:,:,Z(i,:)==1)*T(i));
end


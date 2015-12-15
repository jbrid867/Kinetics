%fit for actin in dextran

data_phi=[0;.0375;.075;.15];
data_F=[30;36.3;41.2;45.2];

k=.38*10^-7;
s=30-(4.4*10^-6)/k;
data_F=k*(data_F-s);

x0=[10^-6;0.95];

x=lsqcurvefit(cfib_actin,x0,data_phi,data_F);


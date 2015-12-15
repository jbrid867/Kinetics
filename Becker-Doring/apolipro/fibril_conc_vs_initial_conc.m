data_c0=[.4;.5;.6;.7];
data_A4=[.0713;.110;.152;.193];
close all


a=10^6;
b=1;

c0=linspace(0.35,0.75,100);
c1_eq=zeros(1,100);
cfib=zeros(1,100);


%scale
data_A4=5.957*data_A4 - 12.46*(data_A4.^2);
scatter(data_c0,data_A4);

lambda=a/b;
for i=1:100
    c1_eq(i)=2*c0(i)/(1+2*lambda*c0(i)+(1+4*lambda*c0(i))^0.5);
    cfib(i)=c0(i)-c1_eq(i)-lambda*c1_eq(i)^2;
end
hold on
plot(c0,cfib);

%fibril % vs initial mass
close all
a=10^4;
b=2*10^-2;
c_0=[7.5;9.6;12.4;14.2;16.2;18.4;20.5].*10^-6;
data=[5;7.6;10.8;12.4;15.3;17.1;19.1].*10^-6;
N=30;

L=length(c_0);
cfib=zeros(L,1);
cr=zeros(N,1);
aob=a/b;

c1=zeros(1,L);
for i=1:L
    c1(i)=2*c_0(i)/(1+2*aob*c_0(i)+(1+4*aob*c_0(i))^0.5);
    cfib(i)=c_0(i)-c1(i);
    for j=2:2
        cfib(i)=cfib(i)-((c1(i)*aob)^(j-1))*c1(i)^2;
    end
   
        
end

plot(c_0,cfib); % plots for fig 5 in wegner actin paper
hold on
scatter(c_0,data);
hold off
% two very small conc data points may not have reached equilibrium

a=10^4;
b=5*10^-2;
c_0=[6.7;8.5;11.5;14.9;17.3;20.3;22.9].*10^-6;
data=[4.7;7.1;9.4;12.6;14.4;16.8;19.4].*10^-6;
N=30;

L=length(c_0);
cfib=zeros(L,1);
cr=zeros(N,1);
aob=a/b;

c1=zeros(1,L);
for i=1:L
    c1(i)=2*c_0(i)/(1+2*aob*c_0(i)+(1+4*aob*c_0(i))^0.5);
    cfib(i)=c_0(i)-c1(i);
    for j=2:2
        cfib(i)=cfib(i)-((c1(i)*aob)^(j-1))*c1(i)^j;
    end
   
        
end
figure
plot(c_0,cfib); % plots for fig 5 in wegner actin paper
hold on
scatter(c_0,data, 's', 'filled');
xlab=xlabel('c_0');
ylab=ylabel('c_{fib} (\muM)');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([tit, xlab, ylab], ...
    'FontName'   , 'AvantGarde');
set([gca]             , ...
    'FontSize'   , 8           );
set([xlab, ylab]  , ...
    'FontSize'   , 10          );
set( tit                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );


hold off

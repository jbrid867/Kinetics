function fib_min_data = becker_fit(x)
% beckdoring script
N=1000;
ICs=zeros(N,1);
ICs(1)=5*10^-6;
ICderiv=zeros(N,1);
phi=0;
% x(1) = disociation constant initial guess (10^-2)
% x(2) = association constant  factor (b*factor=a) initial guess (2.35*10^4)
% x(3) = nucleation guess (3.5*10^4)
% x(4) = data scale initial guess (1.25*10^-7)

data=load('dez_data2');
data(:,1)=data(:,1)*60;

sol = ode15s(@(t,c)beckdoring(t,c,phi,x(1),x(2),x(3)),[0,data(length(data),1)],ICs);

T1=data(:,1);
Y1=deval(sol,T1);
Y1=Y1.';

t=length(T1);
P=zeros(t,1);
M=zeros(t,1);
L=zeros(t,1);
mass=zeros(t,1);
for i=3:N
    P(:)=P(:)+Y1(:,i);
    M(:)=M(:)+i*Y1(:,i);
end
for i=1:N
    mass(:)=mass(:)+i*Y1(:,i);
end
L=M./P;

data(:,2)=(1.28*10^-7)*data(:,2); %scale set to match equilibrium prediction

fib_min_data=zeros(length(data));
fib_min_data=M(:)-data(:,2);


end


function fib_min_data = becker_fit_phi(x)
%function to fit kinetic data with crowders
%   params = [b, factor, k_nuc, scale] not sure if scale is necessary NOT
%   USING CURRENTLY
params=zeros(4,1);
params(1) = 10e-003; %b
params(2) = 2.35000000000000e+006; %factor
params(3) = 35.0000000000000e+003; % k_nuc
params(4) = (1.28*10^-7); %scale???? might need to fit again

%initilal guess for rc is 1.35(nm)

N=1000;
ICs=zeros(N,1);
ICs(1)=5*10^-6;
ICderiv=zeros(N,1);
phi=0.0375;
% params(1) = disociation constant initial guess (10^-2)
% params(2) = association constant  factor (b*factor=a) initial guess (2.35*10^4)
% params(3) = nucleation guess (3.5*10^4)
% params(4) = data scale initial guess (1.25*10^-7)
% x = rc

data=load('act_dex_phi0375.txt');
data(:,1)=data(:,1)*60;

sol = ode15s(@(t,c)beckdoring_rc(t,c,phi,params(1),params(2),params(3),x(1)),[0,(10^4)*data(length(data),1)],ICs); % scaled to meet tolerance

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

data(:,2)=(1.28*10^-7)*data(:,2); %initial scale guess is params(4)

fib_min_data=zeros(length(data));
fib_min_data=M(:)-data(:,2);

end


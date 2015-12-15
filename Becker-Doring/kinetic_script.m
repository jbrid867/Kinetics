%% becker-doring kinetic model. mass conserving

N=100;
ICs=zeros(N,1);
ICs(1)=1;
ICderiv=zeros(N,1);
phi=0;

[T1,Y1] = ode23tb(@(t,c)kinetic_beck(t,c,phi),[0,5],ICs);


t=length(T1);
P=zeros(t,1);
M=zeros(t,1);
L=zeros(t,1);
mass=zeros(t,1);
for i=2:N
    P(:)=P(:)+Y1(:,i);
    M(:)=M(:)+i*Y1(:,i);
end
for i=1:N
    mass(:)=mass(:)+i*Y1(:,i);
end
L=M./P;


plot(T1,M)


phi=0.05;
[T,Y] = ode23tb(@(t,c)kinetic_beck(t,c,phi),[0,5],ICs);

t=length(T);
P=zeros(t,1);
M=zeros(t,1);
L=zeros(t,1);
mass2=zeros(t,1)
for i=2:N
    P(:)=P(:)+Y(:,i);
    M(:)=M(:)+i*Y(:,i);
end
for i=1:N
    mass2(:)=mass2(:)+i*Y(:,i);
end

L=M./P;
hold on
plot(T,M,'r');


phi=0.2;
[T2,Y2] = ode23tb(@(t,c)kinetic_beck(t,c,phi),[0,5],ICs);

t=length(T2);
P=zeros(t,1);
M=zeros(t,1);
L=zeros(t,1);
mass3=zeros(t,1);
for i=2:N
    P(:)=P(:)+Y2(:,i);
    M(:)=M(:)+i*Y2(:,i);
end
for i=1:N
    mass3(:)=mass3(:)+i*Y2(:,i);
end

L=M./P;
plot(T2,M,'c');

figure
hold on
plot(T1,mass);
plot(T,mass2,'r');
plot(T2,mass3,'c');


hold off
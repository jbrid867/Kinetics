function dc = becker(t,c)
N=1000;
dc=zeros(N,1);
J=zeros(N-1,1);
J_1=0;
a=100; % this is the association rate, should be a vector but i will set all to one
b=10^-4;

%build J's
for i=1:N-1
    J(i)=a*c(i)*c(1)-b*c(i+1);
    J_1=J_1+J(i); %sum for c_1
end

dc(1)=(-J(1)-J_1); % special case for c1

for i=2:N-1
    dc(i)=(J(i-1)-J(i)); % main stuff
end
dc(N)=J(N-1); %truncate


end


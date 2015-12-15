function R2 = diffusion( data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pos=zeros(3,1);
pos(1)=data(1);
pos(2)=data(2);
pos(3)=data(3);
R2=zeros(length(data)/3,1);
R2(1)=pos(1)^2+pos(2)^2+pos(3)^2;
for i=2:(length(data)/3)-1
    pos(1)=pos(1)+data(3*i);
    pos(2)=pos(2)+data(3*i+1);
    pos(3)=pos(3)+data(3*i+2);
    R2(i)=pos(1)^2+pos(2)^2+pos(3)^2;

end


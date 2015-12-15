function r = 1oneDdiff( data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
r(1)=data(1);
for i=2:length(data)
    r(i)=r(i-1)+data(i);

end


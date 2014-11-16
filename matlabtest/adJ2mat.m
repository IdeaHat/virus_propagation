function [ A ] = adJ2mat( Al )
%ADJ2MAT Summary of this function goes here
%   Detailed explanation goes here
A = false(numel(Al));

for (i=1:numel(Al))
    A(i,Al{i})=true;
    A(Al{i},i)=true;
end

end


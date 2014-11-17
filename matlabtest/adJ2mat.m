function [ A ] = adJ2mat( Al )
%ADJ2MAT Converts and Adjacency List to an Adjacency Matrix
A = false(numel(Al));

for (i=1:numel(Al))
    A(i,Al{i})=true;
    A(Al{i},i)=true;
end

end


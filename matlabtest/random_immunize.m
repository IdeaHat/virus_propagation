function [ Al ] = random_immunize( Al,k )
%RANDOM_IMMUNIZE Immunizes an adjacency list
%  randomly

immune = randperm(numel(Al),k);
Al(immune)=cell(size(immune));
Al = cellfun(@(x) setdiff(x,immune),Al,'UniformOutput',false);

end


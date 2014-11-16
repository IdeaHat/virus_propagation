function [ Al ] = random_immunize( Al,k )
%RANDOM_IMMUNIZE Summary of this function goes here
%   Detailed explanation goes here
immune = randperm(numel(Al),k);
Al(immune)=cell(size(immune));
Al = cellfun(@(x) setdiff(x,immune),Al,'UniformOutput',false);

end


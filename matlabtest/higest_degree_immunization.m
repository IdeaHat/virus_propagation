function [ Al ] = higest_degree_immunization( Al,k )
%HIGEST_DEGREE_IMMUNIZATION Immunizes an adjacency list
%  based on the highest degree

[~,immune]=sort(cellfun(@(x) numel(x),Al),'descend');
immune = immune(1:k);
Al(immune)=cell(size(immune));
Al = cellfun(@(x) setdiff(x,immune),Al,'UniformOutput',false);

end
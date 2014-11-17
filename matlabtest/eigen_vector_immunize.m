function [ Al ] = eigen_vector_immunize( Al,k,eigen_vec )
%EIGEN_VECTOR_IMMNIZE Immunizes a adjacency list based
%  on the eigen values

if nargin<3
    [~,eigen_vec]=eigs(adj2mat(Al),1);
end

[~,immune]=sort(abs(eigen_vec),'descend');
immune=immune(1:k);
Al(immune)=cell(size(immune));
Al = cellfun(@(x) setdiff(x,immune),Al,'UniformOutput',false);
end


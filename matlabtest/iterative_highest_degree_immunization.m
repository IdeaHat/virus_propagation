function [ Al ] = iterative_highest_degree_immunization( Al,k )
%ITERATIVE_HIGHEST_DEGREE Immunizes an adjacency list
%  based on the highest degree, iteratively removing the
%  highest degree.

A = adJ2mat(Al);
D = cellfun(@(x) numel(x),Al);
while (k>0)
  [~,immune] = max(D);
  
  Al{immune}=[];
  
  Al(A(:,immune)) = cellfun(@(x) setdiff(x,immune),Al(A(:,immune)),'UniformOutput',false);
  D(A(:,immune))=D(A(:,immune))-1;
  D(immune) = 0;
  A(:,immune)=0;
  A(immune,:)=0;
  k = k-1;
end

end


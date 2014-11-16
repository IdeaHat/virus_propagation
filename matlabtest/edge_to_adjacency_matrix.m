function A=edge_to_adjacency_matrix(E,V)
A=zeros(numel(V));
A(sub2ind(size(A),E(:,1),E(:,2)))=1;
A=A|A';

end
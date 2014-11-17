function ADJ=edge_to_adjacency_list(E,V)

ADJ=cell(size(V));
for (i=1:size(E,1))
    ADJ{E(i,1)}=[ADJ{E(i,1)},E(i,2)];
    ADJ{E(i,2)}=[ADJ{E(i,2)},E(i,1)];
end
end
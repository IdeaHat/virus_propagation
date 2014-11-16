function [V,E,graph_size]=read_edge_list(str)
  E = dlmread(str);
  graph_size=E(1,:);
  E=E(2:end,:)+1;
  V=1:graph_size(1);
end
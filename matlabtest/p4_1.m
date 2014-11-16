clear all;
close all;
num_sims=50;
sim_iters=100;
beta = [0.2,0.01];
delta = [0.70,0.60];
k1=200;

file = '~/Downloads/static.network';

[V,E,graph_size] = read_edge_list('~/Downloads/static.network');
Al_base = edge_to_adjacency_list(E,V);
A_base = adJ2mat(Al_base);
base_count = numel(Al_base);
[eig_vec_base,~] = eigs(double(A_base),1);

results = table();

immunization_methods = {...
    @(al,k) al,...
    @(al,k) random_immunize(al,k),...
    @(al,k) higest_degree_immunization(al,k),...
    @(al,k) iterative_highest_degree_immunization(al,k),...
    @(al,k) eigen_vector_immunize(al,k,eig_vec_base)};

immunization_names = {...
    'None',...
    'A: Random', 'B: Highest Degree', 'C: Iterative Highest Degree',...
    'D: Eigen Values'};
for (immval = 1:numel(immunization_methods))
    
    imm_res = table();
    imm_name = immunization_names{immval};
    fprintf('Immunizing with %s...\n',imm_name);
    imm_res.immunization_method = {imm_name};    
    Al = immunization_methods{immval}(Al_base,k1);
    A = adJ2mat(Al);
    
    lambda = eigs(double(A),1);
   
    
    for (i=1:2)
        fprintf('%s param %d...\n',imm_name,i);
        imm_res2 = imm_res;
        imm_res2.parameter_level = i;
        
        imm_res2.effective_strength = calc_sis(A,beta(i),delta(i),lambda);
        imm_res2.epidemic = imm_res2.effective_strength > 1;

        sweep = 0:0.01:1;
        figure;
        plot(sweep,max(lambda)*sweep./delta(i));
        hold on;
        plot([0,1],[1,1],'--');
        title(sprintf('Effective strength vs gamma for delta_%d,\nImmunization: %s',i,imm_name));
        xlabel('gamma (probability)');
        ylabel('Effective Strength');

        imm_re2.beta_threshold=1*delta(i)/max(lambda);

        figure;
        plot(sweep,max(lambda)*beta(i)./sweep);
        hold on;
        plot([0,1],[1,1],'--');
        title(sprintf('Effective strength vs delta for gamma_%d,\nImmunization: %s',i,imm_name))
        xlabel('delta (probability)');
        ylabel('Effective Strength');

        imm_re2.delta_threshold = max(lambda)*beta(i);
        asum = zeros(100,1);
        for (j=1:num_sims)
         asum = asum+simulate_sim(Al,beta(i),delta(i),floor(numel(Al)/10),sim_iters);
        end
        asum = asum./(num_sims*base_count);
        imm_res2.simulation_res = {asum};
        figure;
        plot(asum);
        xlabel('Iteration');
        ylabel('Percent infected');
        title(sprintf('Simulation param_%d\nImmunization: %s',i,imm_name));
        
        results=[results;imm_res2];
    end
end


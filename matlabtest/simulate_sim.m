function [ infected_stats ] = simulate_sim( Al,beta,delta,start_count,steps )
%SIMULATE_SIM Runs a simulation of a SIS infection in a system

n = numel(Al);
infected_stats = zeros(steps,1);
infected_stats(1)=start_count;
infected = false(n,1);
infected(randperm(n,start_count))=true;
infected_count = start_count;
next_infect = false(n,1);
next_cure = false(n,1);
adj = zeros(n,1);
for (iter=2:steps)
    next_infect(:)=false;
    next_cure(:)=false;
    adj(:)=0;
    %find next to infected (O(m))
    for (i=1:n)
        if (infected(i))
            adj(Al{i})=adj(Al{i})+1;
        end
    end
    adj = adj .* ~infected;
    prob = 1-(1-beta).^adj;
    adj_count = sum(adj>0);
    next_infect(adj>0)=rand([adj_count,1])<=beta;
    %next_infect = rand([n,1])<=prob;
    next_cure(infected) = rand([infected_count,1])<=delta;
    infected = (infected & ~next_cure) | next_infect;
    infected_count = sum(infected);
    infected_stats(iter)=infected_count;
end
end
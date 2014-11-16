function [ eff_str ] = calc_sis( A,xmission_prob,healing_prob,max_eigen )
%CALC_SIS Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 4)
    max_eigen = eigs(A,1);
end

eff_str = max_eigen*xmission_prob/healing_prob;

end


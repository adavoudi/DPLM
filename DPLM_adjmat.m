% This function calculates the adjacency matrix (Adj) and 
% the means (M) of every point's neighbourhood (according to the 
% adjacency matrix).
% Note that this is an within class adjacancy matrix (i.e. the
% value of Adj(i,j) for the two points i and j of two different 
% classes is zero.
% 
% 
% Please cite the following paper if you use 
% this code:
%		Davoudi, Alireza, Saeed Shiry Ghidary, and Khadijeh Sadatnejad. 
%		"Dimensionality reduction based on distance preservation to local 
%		mean for symmetric positive definite matrices and its application 
%		in brain–computer interfaces." Journal of Neural Engineering 14.3 
%		(2017): 036019.
% 
%
function [ Adj, M ] = DPLM_adjmat( data, labels, k )

dist = inf * ones(numel(labels));

for i = 1:numel(labels)
    for j = i+1:numel(labels)
        if(labels(i) == labels(j))
            dist(i,j) = distance_riemann(data(:,:,i), data(:,:,j));
            dist(j,i) = dist(i,j);
        end
    end  
end

M = zeros(size(data,1),size(data,1),numel(labels));
Adj = zeros(numel(labels));
for i = 1:numel(labels)
    [~,idx] = sort(dist(i,:), 'ascend');
    Adj(i,idx(1:k)) = 1;
    M(:,:,i) = mean_covariances(data(:,:,idx),'riemann');
end


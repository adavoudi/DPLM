% This is the unsupervised version of the DPLM algorithm which
% doesn't need the labels. Please cite the following paper if you use 
% this code:
%		Davoudi, Alireza, Saeed Shiry Ghidary, and Khadijeh Sadatnejad. 
%		"Dimensionality reduction based on distance preservation to local 
%		mean for symmetric positive definite matrices and its application 
%		in brain–computer interfaces." Journal of Neural Engineering 14.3 
%		(2017): 036019.
% 
% 
% PARAMS:
% 	data: An m*m*n dimensional matrix which contains n SPD matrix of size m*m
% 	dim : The dimensionality of SPD matrices after DM
%	k	: Number of neighbours
%	Adj	: The adjacency matrix calculated by `DPLM_adjmat` function
%	M	: The means matrix calculated by `DPLM_adjmat` function
%
% RETUTNS:
%	U	: The calculated transformation matrix which can be used as 
%		  below to transform an m*m SPD matrix `x` to an dim*dim 
%		  SPD matrix y:
%				y = U'*x*U
%	obj	: Final value of the objective function
%	Adj	: The adjacency matrix calculated by `DPLM_adjmat` function
%	M	: The means matrix calculated by `DPLM_adjmat` function
%
function [U, obj, Adj, M] = uDPLM( data, dim, k, varargin )

if(~isempty(varargin))
    Adj = varargin{1};
    M   = varargin{2};
    if(length(varargin) > 2)
        verbose = varargin{3};
    else
        verbose = 1;
    end
    [U, obj, Adj, M] = DPLM(data, ones(1,size(data,3)), dim, k, Adj, M, verbose);
else
    [U, obj, Adj, M] = DPLM(data, ones(1,size(data,3)), dim, k);
end


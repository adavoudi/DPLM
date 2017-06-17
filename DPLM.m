% This is the supervised version of the DPLM algorithm which
% needs the labels. Please cite the following paper if you use 
% this code:
%		Davoudi, Alireza, Saeed Shiry Ghidary, and Khadijeh Sadatnejad. 
%		"Dimensionality reduction based on distance preservation to local 
%		mean for symmetric positive definite matrices and its application 
%		in brain–computer interfaces." Journal of Neural Engineering 14.3 
%		(2017): 036019.
% 
% PARAMS:
% 	data: An m*m*n dimensional matrix which contains n SPD matrix of size m*m
% 	dim : The dimensionality of SPD matrices after DM
%	labels : A verctor of size n which contains the label of each point
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
function [U, obj, Adj, M] = DPLM( data, labels, dim, k, varargin )

if(~isempty(varargin))
    Adj = varargin{1};
    M   = varargin{2};
    if(length(varargin) > 2)
        verbose = varargin{3};
    else
        verbose = 1;
    end
else
    [ Adj, M ] = DPLM_adjmat( data, labels, k );
    verbose = 1;
end

u0 = randn(size(data, 1), dim);
u0 = orth(u0);

opts.record = verbose;
opts.mxitr  = 1000;
opts.xtol = .1;
opts.gtol = .1;
opts.ftol = .1;
%obj.tau = 1e-3;
%opts.nt = 1;

t = tic;
[U, obj]= OptStiefelGBB(u0, @objfunc, opts);
tsolve = toc(t);

if(verbose)
    [f,~] = objfunc(u0);
    [flast,~] = objfunc(U);
    
    fprintf('Elapsed time: %f\n', tsolve);
    fprintf('Stoped: %s , Init f: %f, Last f: %f \n',obj.msg, f, flast);
end


    function [F, G] = objfunc(u)
        
        F = 0;
        G = 0;
        for i = 1:size(data, 3)
            for j = 1:size(data, 3)
                if(Adj(i,j) == 0)
                    continue;
                end
                
                Ft = distance_ld(data(:,:,j),M(:,:,i))...
                    -distance_ld(u' * data(:,:,j) * u, u' * M(:,:,i) * u);
                F = F + abs(Ft);
                
                G = G -sign(Ft) * (2 * (data(:,:,j) + M(:,:,i)) * u ...
                    * inv( u' * (data(:,:,j) + M(:,:,i)) * u) ...
                    - data(:,:,j) * u * inv(u' * data(:,:,j) * u) ...
                    - M(:,:,i) * u * inv(u' * M(:,:,i) * u));
                
            end
        end
        
    end

end


% This is a supervised dimensionality reduction algorithm proposed by
% Harandi et al. [1].
% 
% PARAMS:
% 	data: An m*m*n dimensional matrix which contains n SPD matrix of size m*m
% 	dim : The dimensionality of SPD matrices after DM
%	labels : A verctor of size n which contains the label of each point 
%	distType: The distance type. Either `riemann` or `ld`
%	Adj	: The adjacency matrix calculated by `Harandi_adjmat` function
%	verbose	: More verbosity if this is 1
%
% RETUTNS:
%	U	: The calculated transformation matrix which can be used as 
%		  below to transform an m*m SPD matrix `x` to an dim*dim 
%		  SPD matrix y:
%				y = U'*x*U
%	obj	: Final value of the objective function
%	Adj	: The adjacency matrix calculated by `Harandi_adjmat` function
%
% References:
% 			[1] Horev, Inbal, Florian Yger, and Masashi Sugiyama. 
%				"Geometry-aware principal component analysis for symmetric positive 
%				definite matrices." Asian Conference on Machine Learning. 2016.
%
function [U, obj, Adj] = Harandi( data, labels, dim, distType, varargin)

verbose = 1;
if(~isempty(varargin))
    Adj = varargin{1};
    if(length(varargin) == 2)
        verbose = varargin{2};
    end
else
    Adj = Harandi_adjmat( labels );
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
            for j = i+1:size(data,3)
                if(Adj(i,j) == 0)
                    continue;
                end
                if(strcmp(distType,'ld'))
                    F = F + Adj(i,j) * ...
                        distance_ld(u' * data(:,:,i) * u, ...
                        u' * data(:,:,j) * u);
                    
                    G = G + Adj(i,j) * ( ...
                        2 * (data(:,:,i) + data(:,:,j)) * u ...
                        * inv( u' * (data(:,:,i) + data(:,:,j)) * u) ...
                        - data(:,:,i) * u * inv(u' * data(:,:,i) * u) ...
                        - data(:,:,j) * u * inv(u' * data(:,:,j) * u));
                    
                elseif(strcmp(distType,'riemann'))
                    F = F + Adj(i,j) * ...
                        distance_riemann(u' * data(:,:,i) * u, ...
                        u' * data(:,:,j) * u);
                    
                    G = G + Adj(i,j) * ( ...
                        4 * (data(:,:,i) * u * inv(u'*data(:,:,i)*u) -...
                        data(:,:,j) * u * inv(u'*data(:,:,j)*u)) *...
                        logm(u'*data(:,:,i)*u*inv(u'*data(:,:,j)*u)));
                end
                
            end
        end
    end

end
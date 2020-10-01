function [ center, asig, distortion ] = cluster_normalized_kmeans( vv, k , normalize)

% function [ asig, distortion ] = cluster_normalized_kmeans( vv, k )
%
% vv is n-by-k matrix whose row i is the embedding of point i in a k-dimensional
%	    space.
%
% if normalize==1 (default value) then first normalizes the rows of vv
% clusters by using kmeans_ortho_multiple with 15 orthogonal and 10 random initial
% centers. This kmeans algorithm provides consistent (stable) results: 
% it returns same clusterings for same input parameters while
% matlab's built-in kmeans function might return different clusterings
%
% uses global options from Spectral Clustering Library (make sure you've run global_options
% before calling this function)

if nargin<3, normalize=1; end

if normalize 
    p = size( vv, 2 );
    vv = vv./repmat( sqrt(sum(vv.*vv,2)), [1, p ]); % project on sphere
end

n_orth = 15;
n_rand = 10;
global KMEANS_THRESHOLD KMEANS_MAX_ITER
threshold=KMEANS_THRESHOLD;
max_iter=KMEANS_MAX_ITER;

[center, asig, distortion] = kmeans_ortho_multiple( vv', k, n_orth, n_rand, threshold, max_iter );


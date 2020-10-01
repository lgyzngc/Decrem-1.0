function [xx, asig_true, k]=gen_2D_circles(numpoints, radius, sigma)

% generates data randomly distributed around circles in 2D
%
% INPUT:
% numpoints - array containing the number of points in each cluster (circle)
%
% radius - array of radiuses of the circles. Also determines the number of
% clusters
%
% sigma: sigma(1) - standard deviation along 1st dimension
%        sigma(2) - standard deviation along 2nd dimension
% Hence, variance is the same for all clusters, but different for two
% dimensions.
%
% OUTPUT
% xx( 2, Ndata ) - data points (Ndata=sum(numpoints))
% asig_true - assignment of points to clusters (1:k), _row_ vector
% k - number of clusters (=length(radius))
%
% EXAMPLE: 
%       xx=gen_2D_circles([100, 200], [5 10], [0.3 0.5]))
%       figure(); subplot(1,2,1); plot(xx(1,:), xx(2,:), 'b.');
%       title('points');
%       sig=2;
%       S=Sexp_from_points(xx, 1/(2*sig^2));
%       subplot(1,2,2); imagesc(S); title('similarity matrix');
%
%
% $ Authors: Marina Meila, Tatiana Maravina
% $ Last Revision: 06-June-2007
% $ Part of SpectraLib_A


k = length(radius);  % number of circles
if length(numpoints) ~= k
     warning('gen_2D_circles:k', 'length(numpoints)<>k=length(radius). Only first k elements of numpoints will be used!')
     numpoints=numpoins(1:k);
end

Ndata=sum(numpoints);
xx = zeros(2,Ndata);
asig_true = zeros(1,Ndata);
start=1;

for iclust = 1:k;
  ncrt=numpoints(iclust);
  [ xxc yyc ] = pol2cart( (1:ncrt)*2*pi/ncrt, radius( iclust ));
  xx(:, start:(start+ncrt-1)) = [ xxc ; yyc ] + randn(2,ncrt).*repmat( [sigma(1); sigma(2)], [1 ncrt]);
  asig_true(:, start:(start+ncrt-1)) = iclust*ones( 1, ncrt );
  start=start+ncrt;
end



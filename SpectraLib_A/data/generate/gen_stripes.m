function [xx, asig_true, k]=gen_stripes (nstripes, width_ratio, npoints, density_ratio, stripeheight)

% function [xx, asig_true, k]=gen_stripes (nstripes, width_ratio, npoints, density_ratio, stripeheight)
% 
% generates 2-by-n matrix xx of stripes pattern with the following properties:
% nstripes - number of dense stripes (same as number of pairs: sparse+dense stripes)
% width_ratio=width_sparse/width_dense   (width_dense is assumed to be 1 unit)
% npoints - # of points per unit of low density stripe
% density_ratio =  [points per unit in dense stripe]/[points per unit in sparse stripe]
% stripeheight (default is 1) - xx(2,:) ranges from 0 to stripeheight ;
%
% EXAMPLE: 
%         [xx, asig_true, k]=gen_stripes (3, 5, 30, 4);
%         % plot results
%         figure(); plot( xx( 1,:), xx( 2, :), '.' ); % just points
%         plot2Dpoints_with_clusters(xx, asig_true, ['.w';'.r';'.g';'.b';'.c';'.m';'.k']);  % mark different
%         clusters with different colors. function from Spectrul Clustering
%         Library
%
%
% $ Authors: Marina Meila, Tatiana Maravina
% $ Last Revision: 06-June-2007
% $ Part of SpectraLib_A

if nargin<5, stripeheight=1; end

sizes = zeros( 1, 2*nstripes );
xx = [];
asig_true = [];
for istripes = 0:nstripes-1;
    sizes( 2*istripes+1 ) = npoints*width_ratio;
    sizes( 2*istripes+2 ) = npoints*density_ratio;
    xorig = (1+width_ratio)*istripes;
    xxtemp = diag( [width_ratio stripeheight ]) * rand( 2, npoints*width_ratio) + repmat( [ xorig; 0], [ 1, npoints*width_ratio ]); 
    [xdum idum ] = sort( xxtemp( 1,:));
    xx = [ xx xxtemp( :, idum ) ];
    xxtemp = diag( [1 stripeheight ]) * rand( 2, density_ratio*npoints)+ repmat( [ xorig+width_ratio; 0], [ 1, density_ratio*npoints ]); 
    [xdum idum ] = sort( xxtemp( 1,:));
    xx = [ xx xxtemp( :, idum ) ];
    asig_true = [ asig_true (2*istripes+1) * ones( 1, npoints*width_ratio) (2*istripes+2)*ones( 1, npoints*density_ratio ) ];
end;

k = 2*nstripes;



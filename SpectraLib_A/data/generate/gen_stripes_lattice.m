function [xx, asig_true, k]=gen_stripes_lattice(nstripes, nxdense, nxsparse, nydense, nysparse, dxdense, dxsparse, dydense)

% generate a stripes pattern with points on a lattice
% 
%
%
% EXAMPLE: 
%         nstripes = 3;
%         nxdense = 30;
%         nxsparse = 90;
%         nydense = 11;
%         nysparse = 4;
%         dxdense = 1;
%         dxsparse = 4;
%         dydense = 1;
%                  
%         [xx, asig_true, k]=gen_stripes_lattice(nstripes, nxdense,nxsparse, nydense, nysparse, dxdense, dxsparse, dydense);
%         % plot results
%         figure(); plot( xx( 1,:), xx( 2, :), '.' ); % just points
%         plot2Dpoints_with_clusters(xx, asig_true, ['.w';'.r';'.g';'.b';'.c';'.m';'.k']);  % mark different
%         clusters with different colors. function from Spectral Clustering
%         Library
%
% $ Authors: Marina Meila, Tatiana Maravina
% $ Part of SpectraLib_A
% $ Last revision: 06-June-2007

dysparse = dydense * (nydense-1)/ (nysparse-1); % match the lowest and highest values for dense and sparse areas

lxdense = nxdense * dxdense + dxsparse-dxdense;
lxsparse = nxsparse * dxsparse;

ndense = nxdense * nydense;
[ xxdense, yydense ] = meshgrid( (1:nxdense)*dxdense, (0:nydense-1)*dydense );
xxd = [ reshape( xxdense, [ 1, ndense ] ); reshape( yydense, [ 1, ndense ] ) ]; 

nsparse = nxsparse * nysparse;
[ xxsparse, yysparse ] = meshgrid( (1:nxsparse)*dxsparse, (0:nysparse-1)*dysparse );
xxsp = [ reshape( xxsparse, [ 1, nsparse ] ); reshape( yysparse, [ 1, nsparse ] ) ]; 

xx = [];
sizes = [];
asig_true = [];

for istripe = 0:nstripes-1;
   xx = [ xx  xxsp + repmat( [istripe*( lxdense + lxsparse ); 0], [ 1, nsparse ] ) xxd + repmat( [istripe*( lxdense + lxsparse ) + lxsparse; 0 ], [ 1, ndense ]) ];
   sizes = [ sizes nsparse ndense ];
   asig_true = [ asig_true repmat( 2*istripe+1, [1, nsparse ]) repmat( 2*istripe+2, [1, ndense ])];
end;

k = 2*nstripes;

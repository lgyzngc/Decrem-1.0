function [component, csizes, clabels] = connected_components(AA)

% AA square matrix, not necessary symmetric
% AA - matrix of size n*n with no connection nodes removed


% Modification of Marina Meila's code
% component at the end is a piecewise constant vector of length n, each 'constant part'
% corresponds to one connected component

[ ii, jj ] = find( tril( AA+AA', -1));
ne = length( ii );
tmp=size(AA);
n=tmp(1);

% find connected components
component = 1:n;

disp( '...begin connected components...' );
for ie = 1:ne
    node2 = ii( ie );
    node1 = jj( ie );
    c1 = GraphComponent( node1, component );
    c2 = GraphComponent( node2, component );
    if c1 ~= c2 
        component( [c1 c2 ] ) = min( c1, c2 );
    end;
end;
disp( '...final path compression...' );

% restore the "pointers"
for ii = 1:n;
    component(ii)=GraphComponent( ii, component );
end;

clabels = unique( component);
ncomp=length(clabels); % number of connected components


csizes = zeros( 1, ncomp ); %sizes of connected components
for ik = 1:ncomp;
    i1 = find( component == clabels(ik) )';
    l1 = length( i1 );
    csizes( ik ) = l1;
end;



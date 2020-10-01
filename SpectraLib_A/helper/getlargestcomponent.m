function ind = getlargestcomponent(S)
% ind - points that form the largest connected component


 [component, csizes, clabels]=connected_components(S); % find isolated points
 % keep only largest component
 [ nmax, lmax ] = max( csizes );
 ind = find (component==clabels(lmax));
 
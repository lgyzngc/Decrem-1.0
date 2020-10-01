function vv = untelescope( vvlarge, sizes )

% function vv = untelescope( vvlarge, sizes )
% 
% compresses each row of vvlarge to a row of sums over the subsets 
% given by sizes
%
% CONSISTENCY: size( vvlarge, 2 ) == sum( sizes )
%
% EXAMPLE: 
% sizes=[2 3 5];
% vvlarge=[1 1 2 2 2 3 3 3 3 3; -1 1 -2 2 2 -3 3 3 -3 3];
% untelescope(vvlarge, sizes)
% ans =
%
%     2     6    15
%     0     2     3

k = length( sizes );

vv = [];
jjold = 0;
for ii = 1:k;
   jj = jjold + sizes( ii );
   vv = [ vv sum( vvlarge( :, jjold+1:jj ), 2 ) ];
   jjold = jj;
end;

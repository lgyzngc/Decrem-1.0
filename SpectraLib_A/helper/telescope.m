function vvlarge = telescope( vv, sizes )

% function vvlarge = telescope( vv, sizes )
% 
% expands vv by duplicating vv(i, :) (i-th row of vv) sizes(i) times for i<=length(sizes)
% rows >length(sizes) are removed
% 
% CONSISTENCY: dim(vv, 1)=length(sizes) (if < - error, if > - just ignores
% extra rows)
%
% EXAMPLE: 
% sizes=[2; 3];
% vv=[1 2 3; 4 5 6];
% telescope(vv, sizes)
% ans =
%
%     1     2     3
%     1     2     3
%     4     5     6
%     4     5     6
%     4     5     6

k = length( sizes );
vvlarge = [];

for ii = 1:k;
   vvlarge = [ vvlarge; repmat( vv( ii, : ), [ sizes( ii ), 1 ] ) ];
end;

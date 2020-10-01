function xx = gen_sym_hypercube( d, unitlength )

% function xx = gen_sym_hypercube( d, unitlength )
% d = nr.ddimensions
% unitlength = size of the cube
% xx( d, 2^d) all the corners of the hypercube in d dimensions, one
%	      point per column
%
%
% $ Authors: Marina Meila
% $ Last Revision: 06-June-2007
% $ Part of SpectraLib_A

n = 2^d;
xx = zeros( d, n );
nn = n/2;

for dd = 1:d;
   ii = repmat( [ ones( 1, nn) zeros( 1, nn ) ], [ 1, 2^(dd-1) ] );
   xx( dd, : ) = ii;
   nn = nn/2;
end;

xx = xx*unitlength;


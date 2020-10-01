function j = GraphComponent( i, component )
%  function j = GraphComponent( i, component )
%
%  returns the component j=1:n of node i in a (disconnected) graph.
%  component( n ) is a vector of "pointers" to the components of nodes
%  in a graph. a component is identified by its lowest number node.
%  component merging is lazy, but a call to GraphComponent restores
%  the pointers before returning the result.

j = component( i );
if( j ~= i );
  component( i ) = GraphComponent( j, component );
end;
j = component( i );


function S=Sexp_from_points (xx, b) 

% computes similarity between two points in R^d based on
% Sxy=exp(-sum(b_i*(x_i - y_j)^2))
%
% INPUT:
% xx - d-by-n matrix with n d-dimensional points
% b  - 1-by-d vector (if a scalar then autumatically duplicated d times)
%
% OUTPUT:
% S - n-by-n similarity matrix
% 
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
% $ Authors: Tatiana Maravina
% $ Last Revision: 06-June-2007
% $ Part of SpectralLib_A


d=size(xx,1);
n=size(xx,2);

if length(b)==1 
    b=repmat(b, [d 1]);
elseif size(b,2)>1 % we need b to be a column vector
    b=b';
end

xxm=xx.*repmat(sqrt(b), [1 n]);
S=squareform(pdist(xxm').^2);
S=exp(-S);

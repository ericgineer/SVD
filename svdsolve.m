% Calculates the SVD of A using the standard textbook method
% NOTE: this function only works on square matrices.
function [U, S, V] = svdsolve(A)
	% Find eigenvalues of A' * A and compute V

	[V, lambda] = eigs(A' * A);


	n = numel(find(lambda != 0)); % Determine number of eigenvalues

	for i = 1:n
		eigvals(i) = lambda(i,i);
	end

	sv = real(sqrt(eigvals));  % singular values are the square root of the eigenvalues

	% Construct singular value matrix

	S = zeros(size(A));

	dim = numel(find(sv > 0));

	D = zeros(dim,dim);
	for i = 1:dim;
		D(i,i) = sv(i);
	end

	S(1:dim,1:dim) = D;

	% Compute U

	dim = min(size(A));

	U = zeros(dim,dim);

	for i=1:dim;
		U(:,i) = 1/sv(i) * A * V(:,i);
	end
end


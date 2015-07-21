% A function to numerically compute the SVD using Jacobi rotations.
% NOTE: This function only works on square matrices.
function [U,S,V]=svd_jacobi(A)
	tolerance = 1.e-8;
	max_iter = 40;

	n=size(A,1);
	U=A;
	V=eye(n);
	for steps=1:max_iter
	  converge = 0;
	  for j=2:n
		for k=1:j-1
		  % compute [alpha gamma;gamma beta]=(k,j) submatrix of U'*U

		  alpha = sum(U(:,k).^2);
		  beta  = sum(U(:,j).^2);
		  gamma = sum(U(:,k) .* U(:,j));
		  
		  converge = max(converge,abs(gamma)/sqrt(alpha*beta));
		  
		  % compute Jacobi rotation that diagonalizes 
		  %    [alpha gamma;gamma beta]
		  if gamma ~= 0
			zeta = (beta-alpha)/(2*gamma);
			t = sign(zeta)/(abs(zeta)+sqrt(1+zeta^2));
		  else
			% if gamma=0, then zeta=infinity and t=0
			t =0;
		  end
		  c = 1 / sqrt(1 + t^2);
		  s = c * t;
		  
		  % update columns k and j of U
		  T = U(:,k);
		  U(:,k) = c*T-s*U(:,j);
		  U(:,j) = s*T+c*U(:,j);
		  
		  % update matrix V of right singular vectors

		  T = V(:,k);
		  V(:,k) = c * T - s * V(:,j);
		  V(:,j) = s * T + c * V(:,j);

		end
	  end
	  if converge < tolerance
		break;
	  end
	end
	if steps >= max_iter
	  fprintf('ERROr: svd_jacobi failed to converge!');
	end

	% the singular values are the norms of the columns of U
	% the left singular vectors are the normalized columns of U
	for j=1:n
	  singvals(j)=norm(U(:,j));
	  U(:,j)=U(:,j)/singvals(j);
	end
	S=diag(singvals);
end

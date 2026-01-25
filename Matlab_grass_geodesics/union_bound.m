function [ub] = union_bound(X, N)

% Inputs
%X -> Constellation or codebook
% N -> Number of rx. antennas

% Output
% ub -> Union bound of the constellation (mind it is valid up to a constant)

[T,M,K] = size(X);
tsum = 0;

for j = 1:K
    for i = 1:j-1
        dp = det(eye(M) - X(:,:,i)'*X(:,:,j) * X(:,:,j)'*X(:,:,i));
        if dp == 0
            dp = 1e-6;
        end
        tsum = tsum + dp^(-N);
    end
end

ub = real(tsum);
end
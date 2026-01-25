function [dM, indx, dmin, indxr, indxc, dsum] = mat_distances(X, distance)

% Inputs
% X -> Codebook or constellation
% d -> Type of distance


% Outputs

% dM -> Matrix with pairwise distance (the diagonal is set to an arbitrary
% number for numeric convenience
% inx -> Index of the minimum pariwise distance [indxr; indxc]
% dmin -> Mimum pairwise distance
% dsum -> Sum of non-diagonal elements in dM


[~,M,K] = size(X);

if nargin == 1
    distance = 'geodesic';
end

dM = zeros(K); 
for i = 1:K
    for j = i+1:K
        [U,S,V] = svd(X(:,:,i)'*X(:,:,j));
        thetak = acos(diag(S));
        if       strcmp(distance,'chordal')
            d = sqrt(sum(sin(thetak).^2));

        elseif   strcmp(distance, 'geodesic')
            d =  sqrt(sum(thetak.^2));

        elseif   strcmp(distance, 'product')
            d = prod(sin(thetak).^2);
        end
        dM(i,j) = d;
        dM(j,i) = d;
    end
end

dsum = sum(dM(:));

dM = real(dM) + 100*eye(K); %Adjust the diagonal for numeric convenience

[~,indx] = min(dM,[],2); % indexes for minimum distances 
                                        % (for each codeword)

[dmin,I] = min(dM(:)); % total minimum distance 

[indxr, indxc] = ind2sub(size(dM),I); % (indxr;indxc) 
                                                      % codeword is the 
                                                      % pair producing the 
                                                      % minimum distance

end
function indx = MLGrass(C,Y)


% Inputs
% C -> Codebook (tensor of size TxMxL)
% Y -> Received signal TxN

% Outputs
% indx -> ML indx


[~,~,L] = size(C);
loglikelihood = zeros(1,L);

for l = 1:L
    Z = C(:,:,l)' * Y;
    valz = sum(abs(Z(:)).^2); %This is the frobenious norm ^2
    loglikelihood(l) = valz;
    % loglikelihood(l) = real(trace(Y'*C(:,:,l)*C(:,:,l)'*Y)); %The other
    % option (slower in matlab)
end

[~,indx] = max(loglikelihood);

end

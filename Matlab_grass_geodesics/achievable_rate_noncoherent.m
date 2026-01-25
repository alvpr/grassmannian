function R_hat = achievable_rate_noncoherent(C, N, rho, numMC)

%ACHIEVABLE_RATE_NONCOHERENT_DIRECT  Direct MC achievable rate 
% Model: Y = X*H + sqrt(M/(T*rho))*W
%   H ~ CN(0,1) (MxN), W ~ CN(0,1) (TxN), codewords equiprobable.
%
% Inputs:
%   C     -> TxMxL codebook or constellation, C(:,:,l)=X_l
%   N     -> number receive antennas
%   rho   -> SNR parameter
%   numMC -> Monte Carlo trials
%
% Output:
%   R_hat : achievable rate in bits/channel-use, R = (1/T) I(X;Y)


[T, M, L] = size(C);
sigma2 = M/(T*rho);
sigma  = sqrt(sigma2);
I_T = eye(T);

acc = 0.0; % accumulates E[ log2 sum_k exp(Lambda_k - Lambda_true) ]

for it = 1:numMC
    % choose transmitted codeword
    ell = randi(L);
    Xell = C(:,:,ell);

    % draw fading and noise
    H = (randn(M,N) + 1j*randn(M,N))/sqrt(2);
    W = (randn(T,N) + 1j*randn(T,N))/sqrt(2);

    % observation
    Y = Xell*H + sigma*W;

    % compute Lambda_k for all k (direct inverses / solves each time)
    Lambda = zeros(L,1);
    for k = 1:L
        Xk = C(:,:,k);
        Rk = Xk*Xk' + sigma2*I_T;         % (T x T)
        invRk = inv(Rk);                 % explicitly invert (as requested)
        quad = real(trace(Y' * invRk * Y));
        Lambda(k) = -quad - N*log(det(Rk));
    end

    % stable log-sum-exp of (Lambda - Lambda_true)
    a = Lambda - Lambda(ell);
    amax = max(a);
    logsum = amax + log(sum(exp(a - amax)));  % natural log
    acc = acc + logsum / log(2);             % convert to log2
end

R_hat = (log2(L) - acc/numMC) / T;
end
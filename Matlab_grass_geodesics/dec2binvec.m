function binVec = dec2binvec(n, K)
    % Validate inputs
    if n < 0 || floor(n) ~= n
        error('First input must be a non-negative integer');
    end
    if K <= 0 || floor(K) ~= K
        error('Second input K must be a positive integer');
    end
    if n > 2^K - 1
        error('Input n is too large to be represented with K bits');
    end

    % Convert to binary string with K bits, MSB-first
    binStr = dec2bin(n, K);            % MSB-first string of length K
    binVec = double(binStr) - '0';     % Convert to numeric vector
end
function codewords = hamming_codewords(ref_codeword, K)
%Return the binary codewords with hamming distance K


N = length(ref_codeword);          % Codeword length
                  % Desired Hamming distance

% Get all combinations of positions to flip (choose K out of N)
positions = nchoosek(1:N, K);      % Each row = positions to flip

% Initialize matrix to hold codewords
num_codewords = size(positions, 1);
codewords = repmat(ref_codeword, num_codewords, 1);

% Flip the bits at the specified positions
for i = 1:num_codewords
    codewords(i, positions(i, :)) = 1;  % Flip 0 to 1 (since reference is all 0s)
end

end
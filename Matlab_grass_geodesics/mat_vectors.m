function [rM] = mat_vectors(D)

%Given a set of vectors (D is a tensor) return the result of computing the riemmaninan
%metric one-to-one

ND = size(D,3);
rM = zeros(ND);

for k = 1:ND
    for m = 1: ND
        rM(k,m) = riem_metric(D(:,:,k), D(:,:,m));
    end
end
    

end
function [dp] = riem_metric(v1, v2)

%The Riemannian product of two vectors v1 and v2 taken in the same
%point U of the Stiefield manifold

dp = trace(v1' * v2);

end
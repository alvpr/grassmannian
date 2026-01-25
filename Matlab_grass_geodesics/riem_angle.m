function [theta] = riem_angle(v1, v2)

%The Riemannian angle of two vectors v1 and v2 with the same
%base point in Stiefield manifold

theta = acos(real(riem_metric(v1,v2)) / (sqrt(riem_metric(v1,v1)) * sqrt(riem_metric(v2,v2))));

end
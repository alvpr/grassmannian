function [d, thetak] = geodesic_distance(Xi, Xj)
%Calculate the Geodesic (Riemannian) distance between Xi and Xj wich are
%points in the Stiefield manifold
% d -> distance
% thetak -> Vector with princpipal angles
[U,S,V] = svd(Xj'*Xi);
thetak = acos(diag(S));
d = sqrt(sum(thetak.^2));

end
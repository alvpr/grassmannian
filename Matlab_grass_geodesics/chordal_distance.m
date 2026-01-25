function [d, thetak] = chordal_distance(Xi, Xj)
%Calculate the chordal distance between Xi and Xj wich are
%points in the Stiefield manifold
% d -> distance
% thetak -> Vector with princpipal angles

[U,S,V] = svd(Xj'*Xi);
thetak = acos(diag(S));
d = sqrt(sum(sin(thetak).^2));

end
function [d] = product_diversity(Xi, Xj)
%Calculate the procut diversity between Xi and Xj wich are
%points in the Stiefield manifold
% d -> product diversity

M = size(Xi,2);
d = det(eye(M) - (Xi')*(Xj)*(Xj')*(Xi));
%this is the product of the squred sines
end
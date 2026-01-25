function [point] = geodesic(U,t,deltaU)
%Calculates the geodesic starting at U with velocity deltaU 
% U is a point in the Stiefield manifold
%delta U is a vector in the horizontal part of the tangent space to the
%Stiefeld manifold at U
%t is the parameter of the curve, for each t a point in the geodic is
%computed. For t = 0 the resulting point is U.

[Q,S,V] = svd(deltaU, "econ");
point = U*V*diag(cos(t*diag(S)))*V' + Q*diag(sin(t*diag(S)))*V';

end
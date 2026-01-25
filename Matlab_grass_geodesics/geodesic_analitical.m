function [point] = geodesic_analitical(U,t,deltaU)
%Calculates the geodesic starting at U with velocity deltaU 
% U is a point in the Stiefield manifold
%delta U is a vector in the horizontal part of the tangent space to the
%Stiefeld manifold at U
%t is the parameter of the curve, for each t a point in the geodic is
%computed. For t = 0 the resulting point is U.

%Here the analitical expression for the conditions given in the paper is
%used (conditions on T=2*M, U and deltaU)
[T,M] = size(U);
D = deltaU(M+1:end,:);%Here is the inferior submatrix
cosM = cos(t/sqrt(M))*eye(M);
sinM = sqrt(M)*D*sin(t/sqrt(M))*eye(M);

point = [cosM; sinM];
end
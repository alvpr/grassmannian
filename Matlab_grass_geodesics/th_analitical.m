function [point] = th_analitical(U,t1,deltaU1,t2,deltaU2)
%Calculates the principal angles using the analitical expression of the
%geodesics for 

% U is a point in the Stiefield manifold
%delta U is a vector in the horizontal part of the tangent space to the
%Stiefeld manifold at U
%t is the parameter of the curve, for each t a point in the geodic is
%1 and 2 denote the indexes respectively

%Here the analitical expression for the conditions given in the paper is
%used (conditions on T=2*M, U and deltaU)
[T,M] = size(U);
sqM = sqrt(M);
D1 = deltaU1(M+1:end,:);%Here is the inferior submatrix
D2 = deltaU2(M+1:end,:);%Here is the inferior submatrix
%This matrix is the hermitic of the first point by the second point
g1hg2 = cos(t1/sqM)*cos(t2/sqM)*eye(M) + M*sin(t1/sqM)*sin(t2/sqM)*D1'*D2;
point = g1hg2;


end
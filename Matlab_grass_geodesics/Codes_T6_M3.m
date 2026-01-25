%% Codes for T=6 M=3 (L = 2,8 saved as Cg2,Cg8)

%% Compute geodesics
T = 6;
M = 3;
maxt = pi/2 * sqrt(M);
cdimen = M*(T-M); % Complex dimensions of the manifold
nbasis = cdimen * 2; %Number of orthogonal vectors in the tangent space of a point

U = [eye(M); zeros(M)]; %inital point

vD = zeros(T, M, nbasis);

%Create the basis
D = weyl_heisenberg_basis(M);

for k = 1:nbasis/2
    vD(:,:,k) = [zeros(M); D(:,:,k)];
    vD(:,:,nbasis/2 + k) = 1j * [zeros(M); D(:,:,k)];
end

for k = 1: nbasis
vD(:,:,k) = vD(:,:,k)./(sqrt(riem_metric(vD(:,:,k),vD(:,:,k))));
end

% Note that the position k and k + nbasis/2 the riem_metric has imaginary
% part

%Now create the list with all the directions using the vectors on the basis
%and their multiplication by -1
vlist = zeros(T, M, 2*nbasis);

for k = 1:nbasis
    vlist(:,:,k) = vD(:,:,k);
    vlist(:,:,nbasis+k) = -1 * vD(:,:,k);
end

% Now with this list compute one geodesic per direction 
ng = 2*nbasis; % Number of geodesics
nt = 10000; %Number of points per geodesic
t = linspace(0,maxt,nt);

G  = zeros(T,M,ng,nt); % G is to store the points resulting from sampling all the geodesics
distances = zeros(ng,nt); %Stores distances from each point in the geodesic to the intial point

for k = 1 : ng
    for n = 1:nt
        G(:,:,k,n) = geodesic(U,t(n),vlist(:,:,k));
        distances(k,n) = geodesic_distance(U, G(:,:,k,n));
    end
end

%Optial, plot the geodesic distance from departing point
% for k = 1:ng
%     figure()
%     plot(t,distances(k,:),".")
%     grid minor
% end



%% Codebook for L = 2 (Case (i) in the paper) just find the point at the end of the first geodesic
L2 = 2;
Cg2 = cat(3, U, G(:,:,1,end));
[dMp, indxp, dminp, indxrp, indxcp, dsump] = mat_distances(Cg2, 'product');
fprintf('Min distance (product) Cg2 dminp = %.4f\n', dminp);
[dMc, indxc, dminc, indxrc, indxcc, dsumc] = mat_distances(Cg2, 'chordal');
fprintf('Min distance (chordal) Cg2 dminc = %.4f\n', dminc);
[dMg, indxg, dming, indxrg, indxcg, dsumg] = mat_distances(Cg2, 'geodesic');
fprintf('Min distance (geodesic) Cg2 dming = %.4f\n', dming);
ub1 = union_bound(Cg2,1);
ub2 = union_bound(Cg2,2);
ub4 = union_bound(Cg2,4);
fprintf('Union bound Cg2 = %.4f, %.4f, %.4f \n\n', ub1, ub2, ub4);

%Bit labelling is trivial


%% Codebook for L = 8 (Case (iv) in the paper) just find the point at the end of the first geodesic

%The indexes are found manually using the script diametral_set

L8 = 8;
imid = floor(nt/2);
Cg8 = zeros(T,M,L8);
group1 = [1,10,19,28];
group2 = [2,11,20,29];
%Here is the code to adjust x
xvals = 0:1:imid-1;
nx = length (xvals);

%%
% This part is to find geodesic sampling, can be skipt is the sampling 
% point is already known
dMpv = zeros(1,nx); %Vector for product distance
dMcv = zeros(1,nx); %Vector for chordal distance
dMgv = zeros(1,nx); %Vector for geodesic distance
ub1 = zeros(1,nx);
ub2 = zeros(1,nx);
ub4 = zeros(1,nx);

for k = 1:nx
    i1 = imid + xvals(k);
    i2 = imid - xvals(k);
    
 
    for l = 1:length(group1)
        Cg8(:,:,l) = G(:,:,group1(l),i1);
        Cg8(:,:,length(group1) + l) = G(:,:,group2(l),i2);
    end

    [~, ~, dMpv(k), ~, ~, ~] = mat_distances(Cg8, 'product');
    [~, ~, dMcv(k), ~, ~, ~] = mat_distances(Cg8, 'chordal');
    [~, ~, dMgv(k), ~, ~, ~] = mat_distances(Cg8, 'geodesic');
    ub1(k) = union_bound(Cg8,1);
    ub2(k) = union_bound(Cg8,2);
    ub4(k) = union_bound(Cg8,4);
end

[ubm1, iubm1] = min(ub1);
[ubm2, iubm2] = min(ub2);
[ubm4, iubm4] = min(ub4);
%%
figure()
plot(xvals, dMpv, 'k-o','MarkerSize', 3,'Linewidth',1.2)
hold on
plot(xvals, dMcv, 'r-o','MarkerSize', 3,'Linewidth',1.2)
hold on
plot(xvals, dMgv, 'b-o','MarkerSize', 3,'Linewidth',1.2)
hold on
% xline(xvals(iubm1), 'r--');
% grid minor
% xline(xvals(iubm2), 'b--');
% grid minor
xline(xvals(iubm4), 'k--');
grid minor

set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
legend('Product diver.', 'Chordal d.', 'Geodesic d.','UB min')
xlabel('Vals for x', 'Interpreter', 'latex');
title('$T=6, M=3, L=8$','Interpreter', 'latex');
legend('$\textrm{DP}\,(\mathcal{X})$', '$\textrm{d}_\textrm{c}\,(\mathcal{X})$', '$\textrm{d}_\textrm{g}\,(\mathcal{X})$','UB$(\mathcal{X}) N=4$ ','Interpreter', 'latex','Location','best')
%%
%With the chosen x select the final codebook

x=1188;
i1 = imid + x;
i2 = imid - x;


for l = 1:length(group1)
    Cg8(:,:,l) = G(:,:,group1(l),i1);
    Cg8(:,:,length(group1) + l) = G(:,:,group2(l),i2);
end

[dMp, indxp, dminp, indxrp, indxcp, dsump] = mat_distances(Cg8, 'product');
fprintf('Min distance (product) Cg8 dminp = %.4f\n', dminp);
[dMc, indxc, dminc, indxrc, indxcc, dsumc] = mat_distances(Cg8, 'chordal');
fprintf('Min distance (chordal) Cg8 dminc = %.4f\n', dminc);
[dMg, indxg, dming, indxrg, indxcg, dsumg] = mat_distances(Cg8, 'geodesic');
fprintf('Min distance (geodesic) Cg8 dming = %.4f\n', dming);
ub1 = union_bound(Cg8,1);
ub2 = union_bound(Cg8,2);
ub4 = union_bound(Cg8,4);
fprintf('Union bound Cg8 = %.4f, %.4f, %.4f \n', ub1, ub2, ub4);

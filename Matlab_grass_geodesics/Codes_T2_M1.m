%% Codes for T=2 M=1 (L = 2,4 saved as Cg2 and Cg4)


%% Compute geodesics
T = 2;
M = 1;
maxt = pi/2 * sqrt(M); %Maximum travel distance of the geodesic
cdimen = M*(T-M); % Complex dimensions of the manifold
nbasis = cdimen * 2; %Number of orthogonal vectors in the tangent space of a point

U = [eye(M); zeros(M)]; %inital point

%This M=1 is a particular case, the vector basis is trival
ut = 1;
D1 = [0;ut];

%Assign the vectors to conform the basis and normalize (so it is
%orthonormal)
vD = zeros (T, M, nbasis);
vD(:,:,1) = D1;

for k = 1:nbasis/2
    vD(:,:,k) = vD(:,:,k)./(sqrt(riem_metric(vD(:,:,k),vD(:,:,k))));
    vD(:,:,nbasis/2 + k) = 1j * vD(:,:,k);
end

%Now create the vector list with all the directions using the vectors on the 
% initial basis and their multiplication by -1
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

%Union bound calculated for antennas N = 1,2,4
ub1 = union_bound(Cg2,1);
ub2 = union_bound(Cg2,2);
ub4 = union_bound(Cg2,4);
fprintf('Union bound Cg2 = %.4f, %.4f, %.4f \n\n', ub1, ub2, ub4);

%Here bit labelling is trivial


%% Codebook for L = 4 (sample the 4 geodesics, this is Case (ii) in the paper)
%Mind that the product diversity is maximum for the geodesic given by D and
%-1*D
%Note that from the product diversity perspective G1 and G2 are closer than
%G3, G4 that's why choosing different indexes for G1 and G2 (but the same
%for G1 and G3
%Here size of Diametral set is D=4

L4 = 4;
imid = floor(nt/2);
Cg4 = zeros(T,M,L4);
%Here is the code to adjust x
xvals = 0:1:imid-1;
nx = length (xvals);


%%
%This part is to find geodesic sampling, can be skipt is the sampling is 
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
    
    Cg4 = cat(3, G(:,:,1,i1), G(:,:,2,i2), G(:,:,3,i1), G(:,:,4,i2));

    [~, ~, dMpv(k), ~, ~, ~] = mat_distances(Cg4, 'product');
    [~, ~, dMcv(k), ~, ~, ~] = mat_distances(Cg4, 'chordal');
    [~, ~, dMgv(k), ~, ~, ~] = mat_distances(Cg4, 'geodesic');
    ub1(k) = union_bound(Cg4,1);
    ub2(k) = union_bound(Cg4,2);
    ub4(k) = union_bound(Cg4,4);
end

[ubm1, iubm1] = min(ub1);
[ubm2, iubm2] = min(ub2);
[ubm4, iubm4] = min(ub4);

figure()
plot(xvals, dMpv, 'k-o','MarkerSize', 3,'Linewidth',1.2)
hold on
plot(xvals, dMcv, 'r-o','MarkerSize', 3,'Linewidth',1.2)
hold on
plot(xvals, dMgv, 'b-o','MarkerSize', 3,'Linewidth',1.2)
hold on
% xline(xvals(iubm1), 'r--');
% grid minor
xline(xvals(iubm2), 'k--');
grid minor
% xline(xvals(iubm4), 'b--');
% grid minor

set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
legend('$\textrm{DP}\,(\mathcal{X})$', '$\textrm{d}_\textrm{c}\,(\mathcal{X})$', '$\textrm{d}_\textrm{g}\,(\mathcal{X})$','UB$(\mathcal{X}) N=2$ ','Interpreter', 'latex','Location','best')
xlabel('Vals for x', 'Interpreter', 'latex');
title('$T=2, M=1, L=4$','Interpreter', 'latex');


%%
%With the chosen x select the final codebook
x = 1959;
i1 = imid + x;
i2 = imid - x;

for l = 1:L4
    Cg4 = cat(3, G(:,:,1,i1), G(:,:,2,i2), G(:,:,3,i1), G(:,:,4,i2));
end

[dMp, indxp, dminp, indxrp, indxcp, dsump] = mat_distances(Cg4, 'product');
fprintf('Min distance (product) Cg4 dminp = %.4f\n', dminp);
[dMc, indxc, dminc, indxrc, indxcc, dsumc] = mat_distances(Cg4, 'chordal');
fprintf('Min distance (chordal) Cg4 dminc = %.4f\n', dminc);
[dMg, indxg, dming, indxrg, indxcg, dsumg] = mat_distances(Cg4, 'geodesic');
fprintf('Min distance (geodesic) Cg4 dming = %.4f\n', dming);
ub1 = union_bound(Cg4,1);
ub2 = union_bound(Cg4,2);
ub4 = union_bound(Cg4,4);
fprintf('Union bound Cg4 = %.4f, %.4f, %.4f \n\n', ub1, ub2, ub4);

%Bit labelling is straighforward
Bg4 = [0,0;0,1;1,1;1,0];
%Here all points are almost equally space, but it bit labelling asign 1 and
%3 00 and 11, respetively




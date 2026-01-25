%% Codes for T=4 M=2 (L = 2,4,8,16 saved as Cg2,Cg4,Cg8,Cg16)

%% Compute geodesics
T = 4;
M = 2;
maxt = pi/2 * sqrt(M); %Maximum travel distance of the geodesic
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

%Normalize
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

%Test geodesic analitical (optional, use the analytical expression of the
%geodesic to compute, only valid for the considered structred
% Gan  = zeros(T,M,ng,nt); % G is to store the points resulting from sampling all the geodesics
% distances_an = zeros(ng,nt); %Stores distances from each point in the geodesic to the intial point
% 
% for k = 1 : ng
%     for n = 1:nt
%         Gan(:,:,k,n) = geodesic_analitical(U,t(n),vlist(:,:,k));
%         distances_an(k,n) = geodesic_distance(U, Gan(:,:,k,n));
%     end
% end

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

%Bit labelling is trivial


%% Codebook for L = 4 (Case (ii) in the paper)
%Use two geodesic with good product diversity and their opposites G1 and
%G2 and G9,G10 and adjust the sampling
%The indexes are found manually using the script diametral_set

L4 = 4;
imid = floor(nt/2);
Cg4 = zeros(T,M,L4);
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
    Cg4= cat(3, G(:,:,1,i1), G(:,:,4,i2), G(:,:,9,i1), G(:,:,12,i2));
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

fs = 18;
lw = 1.5;
ms = 8;

figure()

plot(xvals./imid*pi/4 *sqrt(M), dMpv, 'k-o','MarkerSize', 3,'Linewidth',1.2)
hold on
plot(xvals./imid*pi/4 *sqrt(M), dMcv, 'r-o','MarkerSize', 3,'Linewidth',1.2)
hold on
plot(xvals./imid*pi/4 *sqrt(M), dMgv, 'b-o','MarkerSize', 3,'Linewidth',1.2)
hold on
xline(xvals(iubm2)./imid*pi/4 *sqrt(M), 'k--','Linewidth',1.2);
grid minor


set(gca,'FontSize',fs)
set(gca,'TickLabelInterpreter','latex')
title('$T=4, M=2, L=4$','Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
set(gca,'FontSize',fs)
set(gca,'TickLabelInterpreter','latex')
legend('$\textrm{DP}\,(\mathcal{X})$', '$\textrm{d}_\textrm{c}\,(\mathcal{X})$', '$\textrm{d}_\textrm{g}\,(\mathcal{X})$','UB$(\mathcal{X}) N=2$ ','Interpreter', 'latex','Location','best')

%%
%With the chosen x select the final codebook
x = 1959;
i1 = imid + x;
i2 = imid - x;

for l = 1:L4
    Cg4= cat(3, G(:,:,1,i1), G(:,:,4,i2), G(:,:,9,i1), G(:,:,12,i2));
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
fprintf('Union bound Cg2 = %.4f, %.4f, %.4f \n\n', ub1, ub2, ub4);

% For bit labelling assing 00 and 11 to G1 and G9, but as in the case of
% T2_M1 the distance is almost the same for every point



%% Codebook for L = 8 (Case (iii) in the paper)
%Here L = D
%Use 4 geodesic with good product diversity wich are G1, G4, G6, G7 and
%their opposites
%From product diversity criterion, G1 is closer to all expect G9
%And the same from distance criterion
L8 = 8;
imid = floor(nt/2);
Cg8 = zeros(T,M,L8);

%Here is the code to adjust x
xvals = 0:1:imid-1;
nx = length (xvals);
%%
% This part is to find geodesic sampling, can be skipt is the sampling is 
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
    
    Cg8 = cat(3, G(:,:,1,i1), G(:,:,4,i2), G(:,:,6,i1), G(:,:,7,i2), G(:,:,9,i1), G(:,:,12,i2), G(:,:,14,i1), G(:,:,15,i2));

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
% grid minor
% xline(xvals(iubm4), 'b--');
grid minor

set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
legend('$\textrm{DP}\,(\mathcal{X})$', '$\textrm{d}_\textrm{c}\,(\mathcal{X})$', '$\textrm{d}_\textrm{g}\,(\mathcal{X})$','UB$(\mathcal{X}) N=2$ ','Interpreter', 'latex','Location','best')
xlabel('Vals for x', 'Interpreter', 'latex');
title('$T=4, M=2, L=8$','Interpreter', 'latex');
%%
%With the chosen x select the final codebook
x=0;
i1 = imid + x;
i2 = imid - x;


Cg8 = cat(3, G(:,:,1,i1), G(:,:,4,i2), G(:,:,6,i1), G(:,:,7,i2), G(:,:,9,i1), G(:,:,12,i2), G(:,:,14,i1), G(:,:,15,i2));


[dMp, indxp, dminp, indxrp, indxcp, dsump] = mat_distances(Cg8, 'product');
fprintf('Min distance (product) Cg8 dminp = %.4f\n', dminp);
[dMc, indxc, dminc, indxrc, indxcc, dsumc] = mat_distances(Cg8, 'chordal');
fprintf('Min distance (chordal) Cg8 dminc = %.4f\n', dminc);
[dMg, indxg, dming, indxrg, indxcg, dsumg] = mat_distances(Cg8, 'geodesic');
fprintf('Min distance (geodesic) Cg8 dming = %.4f\n', dming);
ub1 = union_bound(Cg8,1);
ub2 = union_bound(Cg8,2);
ub4 = union_bound(Cg8,4);
fprintf('Union bound Cg8 = %.4f, %.4f, %.4f \n\n', ub1, ub2, ub4);



%% Codebook for L = 16 (Case (iv) in the paper)
%Here use all 16 geodesics
%Same index for  G1, G2, G7, G8 and their opposites
%Different for the other group

L16 = 16;
imid = floor(nt/2);
Cg16 = zeros(T,M,L16);
%Geodesics in the same group can be sampled with the same index wihtouth
%getting 0 in the diversity product
group1 = [1,4,6,7,9,12,14,15];
group2 = [2,3,5,8,10,11,13,16];

%%
% This part is to find geodesic sampling, can be skipt is the sampling is 
% point is already known
xvals = 0:1:imid-1;
nx = length (xvals);

dMpv = zeros(1,nx); %Vector for product distance
dMcv = zeros(1,nx); %Vector for chordal distance
dMgv = zeros(1,nx); %Vector for geodesic distance
ub1 = zeros(1,nx);
ub2 = zeros(1,nx);
ub4 = zeros(1,nx);

for k = 1:nx
    i1 = imid + xvals(k);
    i2 = imid - + xvals(k);
    
    for l = 1:L16
        if ismember(l, group1)
                Cg16(:,:,l) = G(:,:,l,i1);
        else
                Cg16(:,:,l) = G(:,:,l,i2);
        end
    end

    [~, ~, dMpv(k), ~, ~, ~] = mat_distances(Cg16, 'product');
    [~, ~, dMcv(k), ~, ~, ~] = mat_distances(Cg16, 'chordal');
    [~, ~, dMgv(k), ~, ~, ~] = mat_distances(Cg16, 'geodesic');
    ub1(k) = union_bound(Cg16,1);
    ub2(k) = union_bound(Cg16,2);
    ub4(k) = union_bound(Cg16,4);
end

[ubm1, iubm1] = min(ub1);
[ubm2, iubm2] = min(ub2);
[ubm4, iubm4] = min(ub4);

figure()
%Note that here we use the x range that is used in the paper
plot(xvals./imid*pi/4 *sqrt(M), dMpv, 'k-o','MarkerSize', 3,'Linewidth',1.2)
hold on
plot(xvals./imid*pi/4 *sqrt(M), dMcv, 'r-o','MarkerSize', 3,'Linewidth',1.2)
hold on
plot(xvals./imid*pi/4 *sqrt(M), dMgv, 'b-o','MarkerSize', 3,'Linewidth',1.2)
hold on
xline(xvals(iubm2)./imid*pi/4 *sqrt(M), 'k--','Linewidth',1.2);
grid minor

set(gca,'FontSize',fs)
set(gca,'TickLabelInterpreter','latex')
title('$T=4, M=2, L=16$','Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
set(gca,'FontSize',fs)
set(gca,'TickLabelInterpreter','latex')
legend('$\textrm{DP}\,(\mathcal{X})$', '$\textrm{d}_\textrm{c}\,(\mathcal{X})$', '$\textrm{d}_\textrm{g}\,(\mathcal{X})$','UB$(\mathcal{X}) N=2$ ','Interpreter', 'latex','Location','best')



%%
%With the chosen x select the final codebook
x = 1359;
%ubm
% x=1626;
% x=1528; %to optimize ub
% x=1455;
% x = iubm2-1;

i1 = imid + x;
i2 = imid - x;


for l = 1:L16
    if ismember(l, group1)
        Cg16(:,:,l) = G(:,:,l,i1);
    else
        Cg16(:,:,l) = G(:,:,l,i2);
    end
end

[dMp, indxp, dminp, indxrp, indxcp, dsump] = mat_distances(Cg16, 'product');
fprintf('Min distance (product) Cg16 dminp = %.4f\n', dminp);
[dMc, indxc, dminc, indxrc, indxcc, dsumc] = mat_distances(Cg16, 'chordal');
fprintf('Min distance (chordal) Cg16 dminc = %.4f\n', dminc);
[dMg, indxg, dming, indxrg, indxcg, dsumg] = mat_distances(Cg16, 'geodesic');
fprintf('Min distance (geodesic) Cg16 dming = %.4f\n', dming);
ub1 = union_bound(Cg16,1);
ub2 = union_bound(Cg16,2);
ub4 = union_bound(Cg16,4);
fprintf('Union bound Cg16 = %.4f, %.4f, %.4f \n\n', ub1, ub2, ub4);

%To optimize ub
x = 1528;
i1 = imid + x;
i2 = imid - x;
Cg16ub = zeros(T,M,L16);

for l = 1:L16
    if ismember(l, group1)
        Cg16ub(:,:,l) = G(:,:,l,i1);
    else
        Cg16ub(:,:,l) = G(:,:,l,i2);
    end
end


% Bit labeling example for L = 16 points
Bg16 = zeros(L16,log2(L16));
for n = 1:L16/2
    Bg16(n,:) = dec2binvec(n-1,log2(L16));
end

for n = 1:L16/2
    Bg16(L16/2 + n,:) = dec2binvec(L16-n,log2(L16));
end

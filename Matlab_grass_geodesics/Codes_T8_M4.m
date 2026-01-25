%% Codes for T=8 M=4 (L = 2,16,64 saved as Cg2,Cg16,Cg64)

%% Compute geodesics
T = 8;
M = 4;
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

%Codebook for L = 2 (just find the point at the end of the first geodesic
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


%% Codebook for L = 16 (Case (iv) in the paper)

L16 = 16;
imid = floor(nt/2);
Cg16 = zeros(T,M,L16);

%The indexes are found manually using the script diametral_set
group1 = [1,25,16,32, 33,57,48,64];
group2 = [2,20,15,31, 34,52,47,63];

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
    i2 = imid - xvals(k);
    
    for l = 1:length(group1)
        Cg16(:,:,l) = G(:,:,group1(l),i1);
        Cg16(:,:,length(group1) + l) = G(:,:,group2(l),i2);
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
xlabel('Vals for x', 'Interpreter', 'latex');
title('$T=8, M=4, L=16$','Interpreter', 'latex');
legend('$\textrm{DP}\,(\mathcal{X})$', '$\textrm{d}_\textrm{c}\,(\mathcal{X})$', '$\textrm{d}_\textrm{g}\,(\mathcal{X})$','UB$(\mathcal{X}) N=4$ ','Interpreter', 'latex','Location','best')

figure()
plot(xvals, dMpv, 'k-o','MarkerSize', 3,'Linewidth',1.2)

set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
xlabel('Vals for x', 'Interpreter', 'latex');
title('$T=8, M=4, L=16$','Interpreter', 'latex');
legend('$\textrm{DP}$','Interpreter', 'latex','Location','best')
grid minor

%%
%With the chosen x select the final codebook
x = 696;
% x = 444;
% x = iubm1-1;
% x = iubm2-1;
% x = iubm4-1;

i1 = imid + x;
i2 = imid - x;

for l = 1:length(group1)
    Cg16(:,:,l) = G(:,:,group1(l),i1);
    Cg16(:,:,length(group1) + l) = G(:,:,group2(l),i2);
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
fprintf('Union bound Cg16 = %.4f, %.4f, %.4f \n\n', log10(ub1), log10(ub2), log10(ub4));


%% Codebook for L = 64 (Case (v) in the paper)
%Take all geodesics and sample them at the middle

L64 = 64;
Cg64 = zeros(T,M,L64);
imid = floor(nt/2);

for i=1:L64
    Cg64(:,:,i) = G(:,:,i,imid);
end

[dMp, indxp, dminp, indxrp, indxcp, dsump] = mat_distances(Cg64, 'product');
fprintf('Min distance (product) Cg64 dminp = %.8f\n', dminp);
[dMc, indxc, dminc, indxrc, indxcc, dsumc] = mat_distances(Cg64, 'chordal');
fprintf('Min distance (chordal) Cg64 dminc = %.4f\n', dminc);
[dMg, indxg, dming, indxrg, indxcg, dsumg] = mat_distances(Cg64, 'geodesic');
fprintf('Min distance (geodesic) Cg64 dming = %.4f\n', dming);
ub1 = union_bound(Cg64,1);
ub2 = union_bound(Cg64,2);
ub4 = union_bound(Cg64,4);
fprintf('Union bound Cg64 = %.4f, %.4f, %.4f \n\n', log10(ub1), log10(ub2), log10(ub4));


% Bit labelling example
R = log2(64);
Bg64 = zeros(L64,R);

for n = 1:L64/2
    Bg64(n,:) = dec2binvec(n-1,R);
end

for n = 1:L64/2
    Bg64(L64/2 + n,:) = dec2binvec(L64-n,R);
end
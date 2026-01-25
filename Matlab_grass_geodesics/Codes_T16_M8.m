%% Codes for T=16 M=8 (L = 256 saved as Cg)

%% This is a example of Case (v) in the paper, 
% optimize chordal distance for L = 4M^2

%Compute geodesics
T = 16;
M = 8;
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
nt = 1000; %Number of points per geodesic
t = linspace(0,maxt,nt);

G  = zeros(T,M,ng,nt); % G is to store the points resulting from sampling all the geodesics
distances = zeros(ng,nt); %Stores distances from each point in the geodesic to the intial point

for k = 1 : ng
    for n = 1:nt
        G(:,:,k,n) = geodesic(U,t(n),vlist(:,:,k));
        distances(k,n) = geodesic_distance(U, G(:,:,k,n));
    end
end


%%
Lmax = 4*M^2;
Cg = zeros(T,M,Lmax);
imid = floor(nt/2);

for i=1:Lmax
    Cg(:,:,i) = G(:,:,i,imid);
end

[dMp, indxp, dminp, indxrp, indxcp, dsump] = mat_distances(Cg, 'product');
fprintf('Min distance (product) Cg dminp = %.8f\n', dminp);
[dMc, indxc, dminc, indxrc, indxcc, dsumc] = mat_distances(Cg, 'chordal');
fprintf('Min distance (chordal) Cg dminc = %.4f\n', dminc);
[dMg, indxg, dming, indxrg, indxcg, dsumg] = mat_distances(Cg, 'geodesic');
fprintf('Min distance (geodesic) Cg dming = %.4f\n\n', dming);



%% Determine D (size of diamtreal set) and extract the indexes of two of them

%% This is kind of an auxiliar file, it is condidered that G (the tensor of geodesics)
% is already defined (as well as M and T)

%Since the geodesics are already computed, the diametral sets can be
%extracted from them. Looking at the mid points is a straighforward manner
%(and also just the first half of them since the other half is obtaingin
%just by multiplyng by -1 the vector of the geodesic.

%The vectors originating those geodesics are in the same diametral set if
%pairwise the product distance of these points is nonzero

th = 1e-6; % A threshold for numeric pourpuses
Gmid = G(:,:,1:round(end/2),round(end/2));

%Now just obtain the pairwise product diversity matrix
%note this matrix is real and symmetric
K = 2*M^2;
pdmatrix = zeros(K,K);

for k = 1:K
    for m = k+1: K
        res = product_diversity(Gmid(:,:,k), Gmid(:,:,m));
        
        if res < th
            res = 0;
        end

        pdmatrix(m,k) = res;
        pdmatrix(k,m) = res;
    end
end
%%
%Work with pdmatrix to extrac d = D/2 and the indexes

A = (pdmatrix > 0);          % adjacency
A(1:K+1:end) = false;        % ignore diagonal
A = A | A.';                 % enforce symmetry

d = 0;

% Stack holds {cand, sz} where:
% cand = candidates that can still be added
% sz   = current clique size built so far
stack = {1:K, 0};

while ~isempty(stack)
    sz   = stack{end,2};
    cand = stack{end,1};
    stack(end,:) = [];

    % Update best
    if sz > d
        d = sz;
    end

    % Simple bound: even taking all remaining candidates can't beat d
    if sz + numel(cand) <= d
        continue;
    end

    for t = 1:numel(cand)
        v = cand(t);
        next = cand(t+1:end);
        next = next(A(v, next));          % keep only neighbors of v
        stack(end+1,:) = {next, sz+1};    %#ok<SAGROW>
    end
end
%%
%Now find The first diametral set (contains index 1), forcing 2nd index to
%maximize pdmatrix(1, j)

target = d;
best = [];

cand1 = find(A(1,:));
if isempty(cand1)
    best = [];
else
    % pick j that maximizes pdmatrix(1,j) among candidates
    [~,ii] = max(pdmatrix(1,cand1));
    j = cand1(ii);

    % start clique with [1 j]
    clique0 = [1 j];
    cand0 = find(A(1,:) & A(j,:));   % must connect to both 1 and j
    cand0(cand0==1 | cand0==j) = [];

    if target == 1
        best = 1;
    elseif target == 2
        best = clique0;
    else
        stack = {cand0, 2, clique0}; % {cand, sz, clique}

        while ~isempty(stack)
            cand   = stack{end,1};
            sz     = stack{end,2};
            clique = stack{end,3};
            stack(end,:) = [];

            if sz == target
                best = clique;
                break
            end

            if sz + numel(cand) < target
                continue;
            end

            for t = 1:numel(cand)
                v = cand(t);
                next = cand(t+1:end);
                next = next(A(v,next));
                stack(end+1,:) = {next, sz+1, [clique v]}; %#ok<SAGROW>
            end
        end
    end
end

d1indx = best;
%%
%Try to find the second diametral set, disjoint from d1indx, forcing 2nd
%index to maximize pdmatrix(j0, j)

target = d;

avail = true(1,K);
avail(d1indx) = false;

j0 = find(avail, 1, 'first');   % first index not in d1indx

d2indx = [];

if isempty(j0)
    d2indx = [];
else
    cand1 = find(avail & A(j0,:));   % available neighbors of j0

    if isempty(cand1)
        d2indx = [];
    else
        % pick j that maximizes pdmatrix(j0,j) among candidates
        [~,ii] = max(pdmatrix(j0,cand1));
        j = cand1(ii);

        clique0 = [j0 j];
        cand0 = find(avail & A(j0,:) & A(j,:));  % must be avail and connect to both
        cand0(cand0==j0 | cand0==j) = [];

        if target == 1
            d2indx = j0;
        elseif target == 2
            d2indx = clique0;
        else
            stack = {cand0, 2, clique0}; % {cand, sz, clique}

            while ~isempty(stack)
                cand   = stack{end,1};
                sz     = stack{end,2};
                clique = stack{end,3};
                stack(end,:) = [];

                if sz == target
                    d2indx = clique;
                    break
                end

                if sz + numel(cand) < target
                    continue
                end

                for t = 1:numel(cand)
                    v = cand(t);
                    next = cand(t+1:end);
                    next = next(A(v,next));
                    stack(end+1,:) = {next, sz+1, [clique v]}; %#ok<SAGROW>
                end
            end
        end
    end
end

if numel(d2indx) ~= target
    d2indx = [];
    disp("Out of elements to construct the second diametral set")
end

%%
fprintf('Size of diametral set (divided by two is) = %.f\n', d);
d1indx
d2indx

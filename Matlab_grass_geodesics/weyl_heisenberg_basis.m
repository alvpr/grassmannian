function D = weyl_heisenberg_basis(d)
% D = weyl_heisenberg_basis(d)
% Returns D of size d x d x d^2 containing the Weyl–Heisenberg operators.
% General case: D(:,:,k) = X^m Z^n, with m = 0..d-1 (outer), n = 0..d-1 (inner).
% Special case d = 2: ordered as [I, X, Z, XZ] (all real).
%
% NOTE: Not normalized. If you want an orthonormal basis under the
% Hilbert–Schmidt inner product, scale each slice by 1/sqrt(d).

    arguments
        d (1,1) {mustBeInteger, mustBePositive}
    end

    if d == 2
        % --- Special 2x2 (all real) in the order [I, X, Z, XZ]
        I = eye(2);
        X = [0 1; 1 0];
        Z = [1 0; 0 -1];
        XZ = X*Z;
        D = zeros(2,2,4);
        D(:,:,1) = I;
        D(:,:,2) = X;
        D(:,:,3) = Z;
        D(:,:,4) = XZ;
        return
    end

    % --- General d
    omega = exp(2*pi*1i/d);        % d-th root of unity
    X = circshift(eye(d), [0 -1]); % shift operator (j -> j+1 mod d)
    Z = diag(omega.^(0:d-1));      % phase operator

    D = zeros(d, d, d^2);          % auto-promotes to complex as needed

    idx = 1;
    for m = 0:d-1
        Xm = X^m;
        for n = 0:d-1
            Zn = Z^n;
            D(:,:,idx) = Xm * Zn;
            idx = idx + 1;
        end
    end
end
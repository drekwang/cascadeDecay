function rho = Urho0(L_superspace, t, rho0_matrix)
    numEig = sqrt(length(L_superspace));
    U_superspace = expm(1i*L_superspace*t);
    % According to PRB 78, 195410 (2008), we only need the n,m element of
    % each Unm to compute rho(n,m), which means we exclude all exciton
    % mixing terms, which is the same thing as saying no photon emission.
    % Note that this approach doesn't work if there is cross-phasing, i.e.
    % transitions between intermediate exciton states due to spin scatter
    % Compute a new U=exp(-iLt)_{no emission efficiently as below
    U = reshape(diag(U_superspace), numEig, numEig);
    rho = U.*rho0_matrix;
end
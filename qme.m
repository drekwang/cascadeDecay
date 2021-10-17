function drhodt = qme(t, rho, Hd, Lrad, Ldep)
    numEig = length(Hd);
    rho_matrix = reshape(rho, numEig, numEig);
    drhodt_matrix = Hd*rho_matrix - rho_matrix*Hd;
    % Add all the radiative decay Lindblads
    for m = 1:numEig
        for l = m+1:numEig
            L =  squeeze(Lrad(m, l, :, :));
            drhodt_matrix = drhodt_matrix ...
                - 1i/2*((L')*(L)*rho_matrix+rho_matrix*(L')*(L)-2*(L)*rho_matrix*(L'));
        end
    end
    % Add all the dephasing Lindblads
    for l = 1:numEig
        L = squeeze(Ldep(l, l, :, :));
        drhodt_matrix = drhodt_matrix ...
            - 1i/2*((L')*(L)*rho_matrix+rho_matrix*(L')*(L)-2*(L)*rho_matrix*(L'));
    end
    % Divide by i
    drhodt_matrix = drhodt_matrix/1i;
    drhodt = reshape(drhodt_matrix, numEig^2, 1);
end
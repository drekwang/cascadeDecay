function [fidelity, rhoPhotonNorm, minEnergy, classicalFidelity, weightsIdealPaths] = getQuantumFidelity(omegaInput, rInput, ...
    dInput, epsrInput, NelectronsInput, GghybInput, GehybInput, include,...
    numPhotons, numTimesInput, depInput)
    %% This one function does (too) many things.
    % 1) Based on system parameters, compute the level diagram
    % 2) Based on a classical fidelity approach, find the paths with the
    % highest weights
    % 3) Compute the relevant parts of the photon density matrix and spit
    % out the quantum fidelity with a Bell or GHZ state
    
    % Constants
    eVPerHartree = 27.2114;
    nmPerBohr = 0.052918;
    Ang = char(197);
    
    % Hamiltonian for a single emitter. The rows and columns are 00, 01,
    % 10, and 11, where the first and second numbers are the number of
    % electrons in the ground and excited orbital, respectively.
    Hemitter = [0 0 0 0;
                0 omegaInput 0 0;
                0 0 0 0;
                0 0 0 omegaInput];
    % Annhilation operator for the excited orbital
    ae = [0 1 0 0;
        0 0 0 0;
        0 0 0 1;
        0 0 0 0];
    % Annilation operator for the ground orbital
    ag = [0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0];
    % Creation operator for the excited orbital
    aedag = ae';
    % Creation operator for the ground orbital
    agdag = ag';
    r = rInput;
    d = dInput;
    epsr = epsrInput;
    Gghyb = GghybInput;
    Gehyb = GehybInput;


    % Sizes
    Nemitters = size(r,1); % Number of emitters
    Nstates = size(Hemitter, 1); % Number of states per emitter
    % Total dimensionality of the Hilbert space
    dimensions = Nstates^Nemitters; 
    Nelectrons = NelectronsInput;

    % Create annihilation and creation operators in total Hilbert space
    % These are denoted with capital letters (A vs. a)
    Ae = zeros(Nemitters, dimensions, dimensions);
    Ag = Ae;
    Aedag = Ae;
    Agdag = Ae;
    for i = 1:Nemitters
        if i == 1
            Aedummy = ae;
            Agdummy = ag;
            Aedagdummy = aedag;
            Agdagdummy = agdag;
        else
            Aedummy = eye(Nstates);
            Agdummy = eye(Nstates);
            Aedagdummy = eye(Nstates);
            Agdagdummy = eye(Nstates);
        end
        for j = 2:Nemitters
            if i == j
                Aedummy = kron(Aedummy, ae);
                Agdummy = kron(Agdummy, ag);
                Aedagdummy = kron(Aedagdummy, aedag);
                Agdagdummy = kron(Agdagdummy, agdag);
            else
                Aedummy = kron(Aedummy, eye(Nstates));
                Agdummy = kron(Agdummy, eye(Nstates));
                Aedagdummy = kron(Aedagdummy, eye(Nstates));
                Agdagdummy = kron(Agdagdummy, eye(Nstates));
            end
        end
        Ae(i, :, :) = Aedummy;
        Ag(i, :, :) = Agdummy;
        Aedag(i, :, :) = Aedagdummy;
        Agdag(i, :, :) = Agdagdummy;
    end

    % Create normalized dipole transition operator and x,y,z components in the
    % complete Hilbert space
    DipNorm = zeros(Nemitters, dimensions, dimensions);
    DipX = zeros(dimensions, dimensions);
    DipY = zeros(dimensions, dimensions);
    DipZ = zeros(dimensions, dimensions);
    for i = 1:Nemitters
        Ddummy = squeeze(Agdag(i, :, :))*squeeze(Ae(i, :, :)) + squeeze(Aedag(i, :, :))*squeeze(Ag(i, :, :));
        DdummyX = d(i, 1)*(squeeze(Agdag(i, :, :))*squeeze(Ae(i, :, :)) + squeeze(Aedag(i, :, :))*squeeze(Ag(i, :, :)));
        DdummyY = d(i, 2)*(squeeze(Agdag(i, :, :))*squeeze(Ae(i, :, :)) + squeeze(Aedag(i, :, :))*squeeze(Ag(i, :, :)));
        DdummyZ = d(i, 3)*(squeeze(Agdag(i, :, :))*squeeze(Ae(i, :, :)) + squeeze(Aedag(i, :, :))*squeeze(Ag(i, :, :)));
        DipNorm(i, :, :) = Ddummy;
        DipX = DipX + DdummyX;
        DipY = DipY + DdummyY;
        DipZ = DipZ + DdummyZ;
    end

    % Create bare-emitter Hamiltonian in the total Hilbert space
    Hemittertotal = zeros(dimensions, dimensions);
    for i = 1:Nemitters
        if i == 1
            Hdummy = Hemitter;
        else
            Hdummy= eye(Nstates);
        end
        for j = 2:Nemitters
            if i == j
                Hdummy = kron(Hdummy, Hemitter);
            else
                Hdummy = kron(Hdummy, eye(Nstates));
            end
        end
        Hemittertotal = Hemittertotal + Hdummy;
    end

    % Dipole-dipole interaction Hamiltonian in total Hilbert space
    Hdipdiptotal = zeros(dimensions, dimensions);
    for i = 1:Nemitters
        for j = i+1:Nemitters
            scaledipdip = norm(d(i,:))*norm(d(j,:))/(epsr*norm(r(i,:)-r(j,:))^3);
            vectordipdip = dot(d(i,:)/norm(d(i,:)), d(j,:)/norm(d(j,:))) ...
                -3*(dot(d(i,:)/norm(d(i,:)),(r(i,:)-r(j,:))/norm(r(i,:)-r(j,:)))) ...
                *(dot(d(j,:)/norm(d(j,:)),(r(i,:)-r(j,:))/norm(r(i,:)-r(j,:))));
            Jij = scaledipdip * vectordipdip;
            Hdummy = Jij*squeeze(DipNorm(i, :, :))*squeeze(DipNorm(j, :, :));
            Hdipdiptotal = Hdipdiptotal + Hdummy;
        end
    end

    % Hybridization interaction Hamiltonian
    Hhybtotal = zeros(dimensions, dimensions);
    count = 1;
    for i = 1:Nemitters
        for j = i+1:Nemitters
            Hgdummy = Gghyb(1, count)*(squeeze(Agdag(i, :, :))*squeeze(Ag(j, :, :))...
                +squeeze(Agdag(j, :, :))*squeeze(Ag(i, :, :)));
            Hedummy = Gehyb(1, count)*(squeeze(Aedag(i, :, :))*squeeze(Ae(j, :, :))...
                +squeeze(Aedag(j, :, :))*squeeze(Ae(i, :, :)));
            count = count+1;
            Hhybtotal = Hhybtotal + Hgdummy + Hedummy;
        end
    end

    Htotal = Hemittertotal*include(1,1)+Hdipdiptotal*include(1,2)+Hhybtotal*include(1,3);

    % Extract the states that have Nelectrons, since the Hamiltonian conserves
    % the total number of electrons
    % This for-loop generates the labels of each row/col that tells us
    % whether each orbital has 0 or 1 electron. The order follows binary,
    % e.g. for three emitters, the states are 00|00|00, 00|00|01, 00|00|10,
    % ... where each pair corresponds to an emitter, and the first of the
    % pair is the number of electrons in the ground orbital and the second
    % is the number of electrons in the excited orbital.
    for i = 0:dimensions-1
        RowColLabels{i+1} = dec2bin(i);
    end
    count = 1;
    % This for-loop counts how many electrons are in each state
    for i = 0:dimensions-1
        RowColLabel = RowColLabels{i+1};
        numElectrons = 0;
        for j = 1:length(RowColLabel)
            numElectrons = numElectrons + str2num(RowColLabel(j));
        end
        if numElectrons == Nelectrons
            RowColListWithNelectrons(count) = i+1;
            count = count+1;
        end
    end
    % Take the rows/cols that correspond to states with Nelectrons
    HtotalNelectrons = Htotal(RowColListWithNelectrons, RowColListWithNelectrons);
    numEig = length(HtotalNelectrons);
    %     % Uncomment this to set H elements below some threshold to zero
    %     Hmax = max(max(abs(HtotalNelectrons)));
    %     % Set all transition dipole moments below some threshold to zero
    %     for i = 1:length(HtotalNelectrons)
    %         for j = 1:length(HtotalNelectrons)
    %             if abs(HtotalNelectrons(i,j)) < 0.001*Hmax
    %                 HtotalNelectrons(i,j) = 0;
    %             end
    %         end
    %     end
    %     HtotalNelectrons;

    % Diagonalize the Hamiltonian to get the level diagram
    [V, D] = eig(HtotalNelectrons);
    energies = diag(D);
    Hd = D; % Hd is the diagonalized Hamiltonian
    DipXNelectrons = DipX(RowColListWithNelectrons, RowColListWithNelectrons);
    DipXNelectronsEigenbasis = V'*DipXNelectrons*V;
    DipYNelectrons = DipY(RowColListWithNelectrons, RowColListWithNelectrons);
    DipYNelectronsEigenbasis = V'*DipYNelectrons*V;
    DipZNelectrons = DipZ(RowColListWithNelectrons, RowColListWithNelectrons);
    DipZNelectronsEigenbasis = V'*DipZNelectrons*V;
    DipTotalNelectronsEigenbasis = (DipXNelectronsEigenbasis.^2 ...
        +DipYNelectronsEigenbasis.^2+DipZNelectronsEigenbasis.^2).^0.5;
    % Simpler variable name since we'll be using it a lot
    dd = DipTotalNelectronsEigenbasis;

    dmax = max(max(abs(DipTotalNelectronsEigenbasis)));
    % Set all transition dipole moments below some threshold to zero
    for i = 1:length(DipTotalNelectronsEigenbasis)
        for j = 1:length(DipTotalNelectronsEigenbasis)
            if abs(DipTotalNelectronsEigenbasis(i,j)) < 0.01*dmax
                DipTotalNelectronsEigenbasis(i,j) = 0;
            end
        end
    end
    % Check error (print if necessary)
    error = HtotalNelectrons*V - V*D;
    
    %% Now that we have set up the operators and computed the level diagram
    %% we can calculate the entanglement fidelity and efficiency both
    %% classically and quantum mechanically
    % Make Lindbladian operators for radiative decay: 
    % L_{lm}=\sqrt{\gamma^r_{lm}}|m><l|=|d_{lm}||m><l|
    % There is one per transition from higher-energy eigenstate l to
    % lower-energy eigenstate m
    % The output of eig happens to give the eigenenergies in order of
    % increasing energy, so the code below assumes this!
    % First two dimensions cover all l->m transitions
    % Second two dimensions are the actual Lindblad operators
    Lrad = zeros(numEig, numEig, numEig, numEig);
    for m = 1:numEig
        for l = m+1:numEig
            Lrad(m, l, m, l) = dd(m, l);
        end
    end
    % Make Lindbladian operators for dephasing
    % Set dephasing rate to a fraction of the fastest radiative decay rate
    gammad = depInput*dInput(1,1)^2; 
    Ldep = zeros(numEig, numEig, numEig, numEig);
    for l = 1:numEig
        Ldep(l, l, l, l) = sqrt(gammad);
    end
    % Superspace means columns of \rho are stacked to form a vector, so the
    % Liouvillian becomes a matrix
    I = eye(numEig, numEig);
    % Compute Liouvillian in superspace
    L_superspace = kron(I,Hd) - kron(Hd.', I);
    % Add Lindbladian transition operators for radiative decay
    for m = 1:numEig
        for l = m+1:numEig
            J = zeros(numEig, numEig);
            J(m,l) = 1;
            Gamma = dd(m, l)^2;
            L_superspace = L_superspace + 1i*0.5*(2*Gamma*kron(conj(J), J) - Gamma*kron(I, J'*J) ...
                - Gamma*kron(J.'*conj(J), I));
        end
    end
    % Add Lindbladian transition operators for dephasing
    for m = 1:numEig
        J = zeros(numEig, numEig);
        J(m,m) = 1;
        Gamma = gammad;
        L_superspace = L_superspace + 1i*0.5*(2*Gamma*kron(conj(J), J) - Gamma*kron(I, J'*J) ...
            - Gamma*kron(J.'*conj(J), I));
    end
    

    % Propagate density matrix with quantum master equation until
    % \rho_gg = threshold to get the time to integrate to
    rho0_matrix = zeros(numEig, numEig);
    rho0_matrix(numEig, numEig) = 1; % Initialize in doubly excited state
    rho0_superspace = reshape(rho0_matrix, numEig^2, 1);
    tmax = 20;
    tspan = [0 tmax];
    [t, rho] = ode45(@(t, rho) qme(t, rho, Hd, Lrad, Ldep), tspan, rho0_superspace);
    % It'd be nice to set rho_gg_threshold to some constant value, like
    % 0.99, but it turns out in some level structures, not all population
    % goes to the ground state. To account for these cases, we find the
    % asymptote in rho_gg and find the time to get some percent of the way
    % there.
    rho_gg = rho(:,1);
    rho_gg_threshold = 0.98*max(rho_gg);
    thresholdIndex = find(rho_gg > rho_gg_threshold, 1, 'first');
    t = t(1:thresholdIndex);
    
    %% Compute dipole correlation functions->photon density matrix->fidelity
    % Get list of density matrix elements we want to compute (since we don't
    % need to compute all of them to get the fidelity)
    % To get this list, use the classical fidelity code to identify every path
    % and its weight. We'll choose some number of paths with the highest
    % weights, e.g. 2 paths for Bell and GHZ states
    p0Input = diag(rho0_matrix);
    [classicalFidelity, weights, paths, path2Index] = getFidelity(...
        energies, DipTotalNelectronsEigenbasis, p0Input);
    [weightsSorted, weightsSortedIndices] = sort(weights, 'descend');
    [B, I] = sort(weights, 'descend');
    numPhotonsPath1(i) = length(paths{I(1)})-1;
    numPhotonsPath2(i) = length(paths{I(path2Index)})-1;
    weightsIdealPaths = [weights(I(1)) weights(I(path2Index))];
    % Find the minimum difference in energy across the two paths with
    % highest weight & unique decay path
    count = 1;
    for j = [1 path2Index]
        path = paths{I(j)};
        for k = 1:length(path)-1
            photonEnergies(count) = abs(energies(path(k)) - energies(path(k+1)));
            count = count+1;
        end
    end
    count = 1;
    for j = 1:length(photonEnergies)-1
        for k = j+1:length(photonEnergies)
            energyDifferences(count) = abs(photonEnergies(j)-photonEnergies(k));
            count = count+1;
        end
    end
    minEnergy = min(energyDifferences);
    
    numTimes = numTimesInput;
    t = linspace(min(t), max(t), numTimes);
    % Dimensionality of rhoPhoton corresponds to number of superpositioned
    % states in the ideal photon state, e.g. 2 for Bell and GHZ states
    numSuperpositionedStates = 2;
    rhoPhoton = zeros(numSuperpositionedStates, numSuperpositionedStates); 
    % numPhotons corresponds to the number of photons in each superpositioned
    % state in the ideal photon state, e.g. 2 for Bell and 3 for GHZ
    % At the moment, I have not yet figured out how to compute the density
    % matrix for superpositions of photon states with different numbers of
    % photons, e.g. |0> + |11>, i.e. we consider only N-photon states
    idealPaths = {paths{weightsSortedIndices(1)}; paths{weightsSortedIndices(path2Index)}};
    for i = 1:numSuperpositionedStates
        for j = 1:numSuperpositionedStates
            tic
            path1 = idealPaths{i};
            path2 = idealPaths{j};
            % Transition operators: si_dag = \sigma_i^\dag, sj = \sigma_j
            si_dag = zeros(numPhotons, numEig, numEig);
            sj = zeros(numPhotons, numEig, numEig);
            for k = 1:numPhotons
                si_dag(k, path1(k), path1(k+1)) = 1;
                sj(k, path2(k+1), path2(k)) = 1;
            end
            % Solve G_{a,b} for all t for each photon
            % G_{a,b} will have numTime^numPhotons elements
            numTimes = length(t);
            % Hard-coding since for now, we'll only ever deal with Bell or GHZ
            if numPhotons == 2
                G = zeros(numTimes, numTimes);
            elseif numPhotons == 3
                G = zeros(numTimes, numTimes, numTimes);
            else
                fprintf('numPhotons must be 2 or 3 because of lazy hard coding!')
                return
            end
            % Compute the dipole correlation function
            if numPhotons == 2
                for k = 1:numTimes
                    for l = k:numTimes
                        % #15 attempt
                        % Rigorously apply correct order of operations guided
                        % by actual physical intuition and the 2008 paper
                        % It works!
                        t1 = t(k);
                        t2 = t(l);
                        M = Urho0(-L_superspace,t1,rho0_matrix);
                        M = squeeze(sj(1, :, :))*M*squeeze(si_dag(1, :, :));
                        M = Urho0(-L_superspace,t2-t1,M);
                        M = squeeze(sj(2, :, :))*M*squeeze(si_dag(2, :, :));
                        G(k,l) = trace(M);
                    end
                end
            elseif numPhotons == 3
                for k = 1:numTimes
                    for l = k:numTimes
                        for m = l:numTimes
                            t1 = t(k);
                            t2 = t(l);
                            t3 = t(m);
                            M = Urho0(-L_superspace,t1,rho0_matrix);
                            M = squeeze(sj(1, :, :))*M*squeeze(si_dag(1, :, :));
                            M = Urho0(-L_superspace,t2-t1,M);
                            M = squeeze(sj(2, :, :))*M*squeeze(si_dag(2, :, :));
                            M = Urho0(-L_superspace,t3-t2,M);
                            M = squeeze(sj(3, :, :))*M*squeeze(si_dag(3, :, :));
                            G(k,l, m) = trace(M);
                        end
                    end
                end
            end
            % Average  over G_{a,b} to get \rho_{a,b}
            GIntegrate = G;
            for m = 1:numPhotons
                GIntegrate = trapz(t, GIntegrate, numPhotons-m+1);
            end
            rhoPhoton(i,j) = GIntegrate;
            toc
        end
    end
    rhoPhotonNorm = rhoPhoton ./ (rhoPhoton(1,1)+rhoPhoton(2,2));
    ideal = [0.5 0.5; 0.5 0.5];
    fidelity = (trace((rhoPhotonNorm^0.5*ideal*rhoPhotonNorm^0.5)^0.5))^2;
end
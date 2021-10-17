function [classicalFidelity, weights, paths, path2Index] = getFidelity(energies, ...
    DipTotalNelectronsEigenbasis, p0Input)
    %% This script takes a diagonalized Hamiltonian, dipole operator, and
    %% initial state to determine the weights of the irreversible decay paths to
    %% compute the classical fidelity with a desired output photon state (goal)
    numStates = length(energies);
    d= DipTotalNelectronsEigenbasis;
    p0 = p0Input;

    %% Compute all possible paths from an initial state
    % Transition dipole moment threshold for inclusion in path determination
    dth = 0.0;
    % Incrementally add paths starting at the initial state
    paths = {find(p0 > 0)};
    numPaths = 1;
    % 0 means the path is complete, 1 means the path can continue to decay
    incompleteVector = [1]; 
    % While not all paths are complete
    while sum(incompleteVector) > 0
        incompletePathIndices = find(incompleteVector == 1);
        for i = incompletePathIndices
            path = paths{i};
            start = path(end);
            % Find all possible transitions from a given incomplete path
            nextIndices = find(d(start,:) > dth);
            % Next steps must have lower energies
            nextIndices = nextIndices(energies(nextIndices) < energies(start));
            for j = nextIndices
                % If this next step is the first next path found for a
                % given incomplete path, just replace the old path with 
                % this new one with the next step added
                if j == nextIndices(1)
                    paths{i} = [path j];
                    nextNextIndices = find(d(j,:) > dth);
                    nextNextEnergies = energies(nextNextIndices);
                    relativeEnergies = (nextNextEnergies < energies(j));
                    % If there's at least 1 possible next transition, the
                    % path is incomplete
                    if sum(relativeEnergies) > 0
                        incompleteVector(i) = 1;
                    else
                        incompleteVector(i) = 0;
                    end
                % If this next step is NOT the first next path found for a
                % given incomplete path, make a new path
                else
                    numPaths = numPaths + 1;
                    paths{numPaths} = [path j];
                    nextNextIndices = find(d(j,:) > dth);
                    nextNextEnergies = energies(nextNextIndices);
                    relativeEnergies = (nextNextEnergies < energies(j));
                    % If there's at least 1 possible next transition, the
                    % path is incomplete
                    if sum(relativeEnergies) > 0
                        incompleteVector(numPaths) = 1;
                    else
                        incompleteVector(numPaths) = 0;
                    end
                end
            end
        end
    end
    numPaths

    %% Construct matrices of rate constants.
    % Due to Fermi's Golden rule, k_ij \propto d_ij^2, assuming emission into
    % free space with approximately frequency independent density of states.
    % We only consider transitions from higher energy to lower energy
    % kIn (kOut) for a given state correspond to transitions that increase 
    % (decrease) its population, where the given state is lower (higher) in
    % energy than the source (drain)
    energiesRows = ones(numStates, numStates);
    energiesCols = ones(numStates, numStates);
    for i = 1:numStates
        energiesRows(i,:) = energiesRows(i,:).*energies(i);
        energiesCols(:,i) = energiesCols(:,i).*energies(i);
    end
    energiesLessThan = (energiesRows < energiesCols);
    energiesGreaterThan = (energiesRows > energiesCols);
    kIn = d.^2.*energiesLessThan;
    kOut = d.^2.*energiesGreaterThan;

    %% Propagate in time
    tspan = [0 20];
    [t, p] = ode45(@(t, p) odefun(t, p, kIn, kOut), tspan, p0);

    %% Compute the relative flux out from each state to eventually compute the
    % weight of each path
    % Row i, column j is the flux out from state i into state j
    integratedOutFlux = zeros(numStates, numStates);
    integratedOutFluxNorm = zeros(numStates, numStates);
    for i = 1:numStates
        for j = 1:numStates
            integratedOutFlux(i,j) = trapz(t, kOut(i,j)*p(:,i));
        end
        % Normalize each row to get percentages of the flux coming out of a
        % given state to each drain state
        integratedOutFluxNorm(i, :) = integratedOutFlux(i ,:) / ...
            sum(integratedOutFlux(i, :));
    end

    %% Compute the weight of each path
    weights = ones(numPaths, 1);
    for i = 1:numPaths
        path = paths{i};
        for j = 1:length(paths{i})-1
            weights(i, 1) = weights(i, 1) * ...
                integratedOutFluxNorm(path(j), path(j+1));
        end
    end
    coefficients = weights.^0.5;
    coefficientsNorm = coefficients / norm(coefficients);

    %% Compute the classical fideity as the overlap between the output and goal
    % Set weightGoal as 1 for top 2 most weighted paths and check later if they
    % happen to correspond to the desired number of photons on each side
    % and to a unique decay path (no overlapping steps!)
    [weightSorted, indicesSorted] = sort(weights, 'descend');
    coefficientsGoal = zeros(numPaths, 1);
    coefficientsGoal(indicesSorted(1), 1) = 1;
    path1 = paths{indicesSorted(1)};
    for j = 2:length(indicesSorted)
        pathj = paths{indicesSorted(j)};
        % If there's no overlap in the states in the decay path, that's
        % our winner. Break out of the for loop.
        if sum(ismember(path1(2:end-1), pathj(2:end-1))) == 0
            coefficientsGoal(indicesSorted(j), 1) = 1;
            path2Index = j;
            break
        end
    end
    coefficientsGoalNorm = coefficientsGoal / norm(coefficientsGoal);
    classicalFidelity = dot(coefficientsNorm, coefficientsGoalNorm)^2;
end
close all;
clear all;

%% Unit conversions
eVPerHartree = 27.2114;
nmPerBohr = 0.052918;

%% Compute entanglement from a system of N 2-level emitters
% As an example, we have N=3 emitters lined up on the x axis
omegaInput = 0.04; % Energy of the g<->e transition in Hartrees
rInput = [0 0 0; % x, y, z position in Bohr for emitter 1
        1/nmPerBohr 0 0; % emitter 2
        2/nmPerBohr 0 0]; % emitter 3
dInput = [6 0 0; % Transition dipole moment vector in e*Bohr for emitter 1
          6 0 0;
          6 0 0];
epsrInput = 1; % Relative permittivity
% Total number of electrons in the system. We use this number to
% drastically reduce the Hilbert space of the level diagram. This is
% important in the case of many emitters.
NelectronsInput = 3; 
% Hybridization interaction energy between ground state orbitals of
% emitters 1 & 2, 2 & 3, and 1 & 3, respectively. All in Hartree
GghybInput = [0.0 0 0]; 
% Hybridization interaction energy between excited state orbitals
GehybInput = [-0.003, 0, -0.003];
% Set to 1 to include a given sub-Hamiltonian, 0 to not
% In order, the Hamiltonians are the bare emitters, dipole-dipole
% interaction, and hybridization interaction
include = [1 1 1];
% The number of photons in each state in the density matrix
% The code is currently hard-coded to take numPhotons = 2 or 3 since in the
% paper we only study Bell states and 3-photon GHZ states, although the
% code can probably be straightforwardly generalized.
numPhotons = 3;
% To compute elements of the density matrix, we integrate over a tensor 
% numPhotons dimensions and numTimesInput elements per dimension
% Increasing numTimesInput gives more converged results, but the
% computation time scales quite poorly.
numTimesInput = 25;
depInput = 0;

% 'fidelity' is the quantum fidelity with a 3-photon GHZ state for
% numPhotons = 3 and Bell state for numPhotons = 2. 
% 'rhoPhotonNorm' is the 2x2 normalized density matrix
% 'minEnergy' is the minimum energy between photons in the desired state
% 'classicalFidelity' is the fidelity computed using the classical method 
% of rate equations that ignores off-diagonal decoherences
% 'weightsIdealPaths' are the populations associated with the decay paths
% involved in the photon state of interested. 
% Summing them gives the efficiency eta in the paper.
[fidelity, rhoPhotonNorm, minEnergy, classicalFidelity, weightsIdealPaths] = ...
    getQuantumFidelity(omegaInput, rInput, dInput, epsrInput, ...
   NelectronsInput, GghybInput, GehybInput, include, numPhotons, numTimesInput, depInput)

%% Below is a script for looping through some array of values to generate
%% plots like those in Fig 2 and 3 of the paper
% Default values. All in atomic units.
omegaDefault = 1/eVPerHartree;
rDefault = 1/nmPerBohr;
dDefault = 6;
epsrDefault = 1;
NelectronsDefault = 3;
GghybDefault = 0.0;
GehybDefault = -0.0002;
includeDefault = [1 1 1];
numPhotonsDefault = 3;
numTimesInput = 25;
depDefault = 0;

% As an example, let's sweep the relative dephasing rate
sweepArray = linspace(0, 2, 2);
rhoPhotonNormArray = zeros(length(sweepArray),2,2);
weightsIdealPaths = zeros(length(sweepArray),2);
fidelityArray = zeros(length(sweepArray), 1);
minEnergyDifferences = zeros(length(sweepArray), 1);
classicalFidelityArray = zeros(length(sweepArray), 1);
for i = 1:length(sweepArray)
    omegaInput = omegaDefault;
    rInput = [0 0 0;
            rDefault 0 0;
            2*rDefault 0 0];
    dInput = [dDefault 0 0;
              dDefault*cos(0.975) dDefault*sin(0.975) 0;
              dDefault 0 0];
    epsrInput = epsrDefault;
    NelectronsInput = NelectronsDefault;
    GghybInput = [GghybDefault 0 GghybDefault];
    GehybInput = [GehybDefault 0 GehybDefault];
    include = includeDefault;
    numPhotons = numPhotonsDefault;
    depInput = sweepArray(1,i);
    [fidelityArray(i), rhoPhotonNormArray(i,:,:), ...
        minEnergyDifferences(i), classicalFidelityArray(i), ...
        weightsIdealPaths(i,:)] = getQuantumFidelity(omegaInput, rInput, ...
    dInput, epsrInput, NelectronsInput, GghybInput, GehybInput, include,...
    numPhotons, numTimesInput, depInput);
end

% Plot
figure()
yyaxis left
plot(sweepArray, abs(fidelityArray), 'Linewidth', 2)
% Quantum fidelity \mathcal{F}
hold on
% Efficiency eta
plot(sweepArray, abs(sum(weightsIdealPaths,2)), '-.', 'Linewidth', 2)
set(gca,'XMinorTick','on','YMinorTick','on')
ylim([0.9*min([abs(min(fidelityArray)) abs(min(sum(weightsIdealPaths,2)))]) ...
    1.2*max([abs(min(fidelityArray)) abs(max(sum(weightsIdealPaths,2)))])])
xlabel('\gamma_d/\gamma_0')
ylabel('$\mathcal{F}$, $\eta$', 'Interpreter', 'latex')
yyaxis right
% Minimum energy difference \Delta E_{min}
plot(sweepArray, minEnergyDifferences*eVPerHartree*1000, '-o','Linewidth', 2)
xlabel('\gamma_d/\gamma_0')
ylabel('\Delta E_{min} [meV]')
set(gca,'XMinorTick','on','YMinorTick','on')
formatFigure()
% Export the data
set(gcf, 'PaperPosition', [0 0 8 16])    % can be bigger than screen 
set(gcf, 'PaperSize', [8 16]);    % Same, but for PDF output
print(gcf, 'fig_GHZ_dephasing.pdf', '-dpdf', '-r300' );

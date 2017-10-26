function arkin98_lambda_bifurcation_ssa
% Lambda phage-infected E. coli lysis vs lysogeny bifurcation stochastic model
%   Based on Arkin, Ross, and McAdams (1998) Stochastic kinetic analysis of
%   developmental pathway bifurcation in phage \lambda-infected Escherichia coli
%   cells. Genetics 149:1633-48
% The original paper (and the papers it references) aren't enough to perfectly
%   reproduce their model so I've filled in the gaps
% The model is pretty complex, especially the transcription and translation
%   "reactions" which include sequential elongation (stepping along nucleotide)
%   steps. The states for major species like proteins are represented using the
%   usual matrix of states x times, but the individual positions of RNAP on DNA
%   and ribosome on mRNA are represented in a sparse method that isn't saved.

clear; close all; clc
rng('default');

% Units used here
%   Time: s
%   Volume: L
%   Conc: M (molar), assuming molar for m insteal of molal or whatever
% In SSA model, convert these to the equivalent per-count values using nA and V

%% Rate constants
k0 = 4.76e-18; % (L)
k1 = 0.0007; % (1/s)
k2 = 0.05; % (1/(M*s))
k3 = 0.5; % (1/s)
k4 = 0.0025; % (1/s)
k5 = 0.05; % (1/(M*s))
k6 = 0.5; % (1/s)
k7 = 0.00231; % (1/s)
k8 = 0.01; % (1/(M*s))
k9 = 0.01; % (1/s)
k10 = 0.002; % (1/s)
k11 = 0.01; % (1/(M*s))
k12 = 0.001; % (1/s)
k13 = 0.0001; % (1/s)
k14 = 0.00025; % (1/(M*s))
k15 = 0.065; % (1/s)
k16 = 0.6; % (1/s)
k17 = 0.01; % (1/(M*s))
k18 = 0.01; % (1/s)
k19 = 0.001; % (1/s)

% Initial compartment volume
V = 1e-15; % (L)

% Constants
nA = 6.022e23; % (count/mol) Avogadro's number
R = 1.9872e-3; % (kcal/(K*mol)) gas constant
T = 273 + 37; % (K) temperature

%% Promoter state energy weights
% P_RE
g = zeros(4,1); % energies
g(1) = 0.0; % reference energy
g(2) = -9.9; % (kcal/mol)
g(3) = -9.7;
g(4) = -21.5;
pg_P_RE = exp(-g./(R*T));

% P_L
g = zeros(10,1); % energies
g(1) = 0.0; % reference energy
g(2) = -10.9; % (kcal/mol)
g(3) = -12.1;
g(4) = -11.7;
g(5) = -10.1;
g(6) = -12.5;
g(7) = -22.9;
g(8) = -20.9;
g(9) = -22.8;
g(10) = 23.7;
pg_P_L = exp(-g./(R*T));

% P_R and P_RM
dG1  = -11.7; % (kcal/mol)
dG2  = -10.1;
dG3  = -10.1;
dG1p = -10.8;
dG2p = -10.8;
dG3p = -12.1;
dGRM = -11.5;
dGR  = -12.5;
dG12 = -1.9; % calculated
dG23 = -2.0; % calculated

g = zeros(40,1); % energies
g(1) = 0.0; % reference energy
g(2) = dG1; % (kcal/mol)
g(3) = dG2;
g(4) = dG3;
g(5) = dG1p;
g(6) = dG2p;
g(7) = dG3p;
g(8) = dGRM;
g(9) = dGR;
g(10) = dG1 + dG2 + dG12;
g(11) = dG1 + dG3;
g(12) = dG2 + dG3 + dG23;
g(13) = dG1p + dG2p;
g(14) = dG1p + dG3p;
g(15) = dG2p + dG3p;
g(16) = dGR + dGRM;
g(17) = dG1 + dG2p;
g(18) = dG1p + dG2;
g(19) = dG1p + dG3;
g(20) = dG1 + dG3p;
g(21) = dG2p + dG3;
g(22) = dG2 + dG3p;
g(23) = dGR + dG3;
g(24) = dG2 + dGRM;
g(25) = dG1 + dGRM;
g(26) = dGR + dG3p;
g(27) = dG2p + dGRM;
g(28) = dG2p + dGRM;
g(29) = dG1 + dG2 + dG3 + dG12;
g(30) = dG1p + dG2p + dG3p;
g(31) = dG1 + dG2 + dG3p + dG12;
g(32) = dG1 + dG2p + dG3;
g(33) = dG1p + dG2 + dG3 + dG23;
g(34) = dG1p + dG2p + dG3;
g(35) = dG1p + dG2 + dG3p;
g(36) = dG1 + dG2p +dG3p;
g(37) = dG1 + dG2 + dGRM + dG12;
g(38) = dG1p + dG2p + dGRM;
g(39) = dG1 + dG2p + dGRM;
g(40) = dG1p + dG2 + dGRM;
pg_P_RRM = exp(-g./(R*T));

%% Promoter initializables rates
P_RE_init_inds = [2,4];
P_RE_init_rates = [0.00004,0.015]; % (1/s)

P_L_init_inds = [6];
P_L_init_rates = [0.011]; % (1/s)

P_RRM_init_inds = [8:9,23:28,37:40];
P_RM_init_rate = 0.011; % (1/s) Shea 85, Table 3, Stimulated P_RM k_Promoter, single val
P_R_init_rate =0.014; % (1/s) Shea 85, Table 3, Basal P_R k_Promoter treated as stimulated, single val

%% Initial conditions - counts
% Main states: x = [Cro, CI, CII, CIII, P1, P2, N, Cro2, CI2, P1_CII, P1_CIII, P2_CII, P2_CIII, RNAP, ribo]
Cro_0  = 0;
CI_0   = 0;
CII_0  = 0;
CIII_0 = 0;
P1_0   = 35e-9; % (M)
P2_0   = 140e-9; % (M)
N_0    = 0;
Cro2_0 = 0;
CI2_0  = 0;
P1_CII_0  = 0;
P1_CIII_0 = 0;
P2_CII_0  = 0;
P2_CIII_0 = 0;
RNAP_0 = 30e-9; % (M) free RNA polymerase
ribo_0 = 500e-9; % (M) free ribosomes

c0 = [Cro_0, CI_0, CII_0, CIII_0, P1_0, P2_0, N_0, Cro2_0, CI2_0, P1_CII_0, P1_CIII_0, P2_CII_0, P2_CIII_0, RNAP_0, ribo_0];

x0 = round(c0*nA*V);

%% Set simulation conditions
tf = 35*60; % (s)

% Run stochastic simulation
%   Each run implicitly uses a different RNG seed
nRuns = 5;
ts = cell(nRuns,1);
xs = cell(nRuns,1);
for i = 1:nRuns
    [tSsa, xSsa] = ssa_sim(tf, x0);
    ts{i} = tSsa';
    xs{i} = xSsa';
end

% Plot trajectories
figure
hold on
ax = gca;
for i = 1:nRuns
    ax.ColorOrderIndex = 1;
    stairs(ts{i}/60, xs{i}) % using regular plot here is not as accurate
end
hold off
xlabel('Time')
ylabel('Counts')
title(sprintf('Arkin 98 lambda Bifurcation'))
xlim([0, tf/60]) % min
% legend('A','B','C','Location','best')


    function [ts, xs] = ssa_sim(tf, x0)
        % Gillespie stochastic simulation algorithm
        % States: x1 = A, x2 = B, x3 = C
        % Rxns: a1 = on, a2 = off
        nt = 1e3; % preallocate times to start with
        ts = zeros(1,nt);
        it = 1;
        t = 0;
        ts(it) = t;
        
        nx = length(x0);
        xs = zeros(nx,nt); % preallocate solution
        x = x0;
        xs(:,it) = x;
        
        while t <= tf
            % Increment solution index
            it = it + 1;
            
            % Reallocate more solution space if needed (double each time)
            if it == nt
                ts = [ts, zeros(1,nt)];
                xs = [xs, zeros(nx,nt)];
                nt = nt * 2;
            end
            
            % Split fixed states for easier use
            Cro  = x(1);
            CI   = x(2);
            CII  = x(3);
            CIII = x(4);
            P1   = x(5);
            P2   = x(6);
            N    = x(7);
            Cro2 = x(8);
            CI2  = x(9);
            P1_CII  = x(10);
            P1_CIII = x(11);
            P2_CII  = x(12);
            P2_CIII = x(13);
            RNAP = x(14);
            ribo = x(15);
            
            % Calculate volume
            %   Cell grows continuously from t = 0 to t = 35 min, doubling in size
            V = (1 + k0*t) * 1e-15; % (L)
            
            %% Calculate reaction propensities
            a = zeros(1,23); % total fixed reactions, additional ones based on txtl will be appended as needed
            a(1) = k1*CI;
            a(2) = k2/(nA*V)*CI*(CI-1); % check /2
            a(3) = k3*CI2;
            a(4) = k4*Cro;
            a(5) = k5/(nA*V)*Cro*(Cro-1); % check /2
            a(6) = k6*Cro2;
            a(7) = k7*N;
            a(8) = k8/(nA*V)*CII*P1;
            a(9) = k9*P1_CII;
            a(10) = k10*P1_CII;
            a(11) = k11/(nA*V)*CIII*P1;
            a(12) = k12*P1_CIII;
            a(13) = k13*P1_CIII;
            a(14) = k14/(nA*V)*CII*P2;
            a(15) = k15*P2_CII;
            a(16) = k16*P2_CII;
            a(17) = k17/(nA*V)*CIII*P2;
            a(18) = k18*P2_CIII;
            a(19) = k19*P2_CIII;
                    
            %% Calculate promoter state probabilities, select, and add rxns
            % Fast promoter state equilibrium assumption
            % If an initiatiable promoter state is selected, add that initiation
            %   to the list of possible reactions (the rest are assumed 0)
            % Initially, when there's no proteins, RNAP binding is heavily favored
            % But the double RNAP P_RM+P_R can't initiate
            P_RE_probs = get_P_RE_probs()';
            P_L_probs = get_P_L_probs()';
            P_RRM_probs = get_P_RRM_probs()';
            
            P_RE_ind = pick_rxn(P_RE_probs);
            P_L_ind = pick_rxn(P_L_probs);
            P_RRM_ind = pick_rxn(P_RRM_probs);
            
            aP = zeros(4,1); % rxn propensities for the 4 promoter sites in order [RE,L,RM,R]
            P_RE_mask = ismember(P_RE_init_inds, P_RE_ind);
            if any(P_RE_mask)
                aP(1) = P_RE_init_rates(P_RE_mask)*RNAP; % kOC is 1st order - is kOC*RNAP the right rate form?
            end
            P_L_mask = ismember(P_L_init_inds, P_L_ind);
            if any(P_L_mask)
                aP(2) = P_L_init_rates(P_L_mask)*RNAP;
            end
            P_RRM_mask = ismember(P_RRM_init_inds, P_RRM_ind);
            if any(P_RRM_mask)
                aP(3) = P_RM_init_rate*RNAP; % both possible directions
                aP(4) = P_R_init_rate*RNAP;
            end
            
            a(20:23) = aP; % add on promoter state rxns
            
            %% Reaction propensity normaliation
            ao = sum(a);
            
            % Calculate next reaction time
            tau = (1/ao)*log(1/rand());
            
            % Update time
            t = t + tau;
            
            % Pick reaction - rxn reaches equilibrium so always picks u = 1 or 2
            u = pick_rxn(a/ao);
            
            % Handle corner case of no rxns possible (i.e., depletion)
            if isinf(t) || isempty(u)
                ts(it) = tf;
                xs(:,it) = x;
                break
            end
            
            % Update species
            % If it's not one of these, then go to sparse rxns
            switch u
                case 1 % CI -> 0
                    CI = CI - 1;
                case 2 % 2CI -> CI2
                    CI = CI - 2;
                    CI2 = CI2 + 1;
                case 3 % CI2 -> 2CI
                    CI = CI + 2;
                    CI2 = CI2 - 1;
                case 4 % Cro -> 0
                    Cro = Cro - 1;
                case 5 % 2Cro -> Cro2
                    Cro = Cro - 2;
                    Cro2 = Cro2 + 1;
                case 6 % Cro2 -> 2Cro
                    Cro = Cro + 2;
                    Cro2 = Cro2 - 1;
                case 7 % N -> 0
                    N = N - 1;
                case 8 % CII + P1 -> P1_CII
                    CII = CII - 1;
                    P1 = P1 - 1;
                    P1_CII = P1_CII + 1;
                case 9 % P1_CII -> CII + P1
                    CII = CII + 1;
                    P1 = P1 + 1;
                    P1_CII = P1_CII - 1;
                case 10 % P1_CII -> P1
                    P1_CII = P1_CII - 1;
                    P1 = P1 + 1;
                case 11 % CIII + P1 -> P1_CIII
                    CIII = CIII - 1;
                    P1 = P1 - 1;
                    P1_CIII = P1_CIII + 1;
                case 12 % P1_CIII -> CIII + P1
                    CIII = CIII + 1;
                    P1 = P1 + 1;
                    P1_CIII = P1_CIII - 1;
                case 13 % P1_CIII -> P1
                    P1_CIII = P1_CIII - 1;
                    P1 = P1 + 1;
                case 14 % CII + P2 -> P2_CII
                    CII = CII - 1;
                    P2 = P2 - 1;
                    P2_CII = P2_CII + 1;
                case 15 % P2_CII -> CII + P2
                    CII = CII + 1;
                    P2 = P2 + 1;
                    P2_CII = P2_CII - 1;
                case 16 % P2_CII -> P2
                    P2_CII = P2_CII - 1;
                    P2 = P2 + 1;
                case 17 % CIII + P2 -> P2_CIII
                    CIII = CIII - 1;
                    P2 = P2 - 1;
                    P2_CIII = P2_CIII + 1;
                case 18 % P2_CIII -> CIII + P2
                    CIII = CIII + 1;
                    P2 = P2 + 1;
                    P2_CIII = P2_CIII - 1;
                case 19 % P2_CIII -> P2
                    P2_CIII = P2_CIII - 1;
                    P2 = P2 + 1;
                case 20 % P_RE initiation
                    
                case 21 % P_L initiation
                    
                case 22 % P_RM initiation
                    
                case 23 % P_R initiation
                    
            end
            
            % Recombine fixed states
            x = [Cro, CI, CII, CIII, P1, P2, N, Cro2, CI2, P1_CII, P1_CIII, P2_CII, P2_CIII]';
            
            % Store results
            ts(it) = t;
            xs(:,it) = x;
        end
        
        % Trim solution
        ts = ts(1:it);
        xs = xs(:,1:it);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helper functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function probs = get_P_RE_probs()
            % Calculate P_RE promoter state probabilities
            %   See Arkin 98 Table 1
            c = zeros(4,1); % concentrations - counts are fine since conversion to concs is just a constant that's factored out
            c(1) = 1; % empty, no concs
            c(2) = RNAP; % (counts)
            c(3) = CII;
            c(4) = CII*RNAP;
            
            weights = pg_P_RE.*c;
            Z = sum(weights);
            probs = weights / Z;
        end
        
        function probs = get_P_L_probs()
            % Calculate P_L promoter state probabilities
            %   See Arkin 98 Table 1
            c = zeros(10,1);
            c(1) = 1; % empty, no concs
            c(2) = Cro2; % (counts)
            c(3) = Cro2;
            c(4) = CI2;
            c(5) = CI2;
            c(6) = RNAP;
            c(7) = Cro2^2;
            c(8) = Cro2*CI2;
            c(9) = CI2*Cro2;
            c(10) = CI2^2;
            
            weights = pg_P_L.*c;
            Z = sum(weights);
            probs = weights / Z;
        end
        
        function probs = get_P_RRM_probs()
            % Calculate P_R and P_RM promoter state probabilities
            %   See Shea 85 Table 2
            c = zeros(40,1);
            c(1) = 1; % empty, no concs
            c(2) = CI2; % (counts)
            c(3) = CI2;
            c(4) = CI2;
            c(5) = Cro2;
            c(6) = Cro2;
            c(7) = Cro2;
            c(8) = RNAP;
            c(9) = RNAP;
            c(10) = CI2^2;
            c(11) = CI2^2;
            c(12) = CI2^2;
            c(13) = Cro2^2;
            c(14) = Cro2^2;
            c(15) = Cro2^2;
            c(16) = RNAP^2;
            c(17) = CI2*Cro2;
            c(18) = CI2*Cro2;
            c(19) = CI2*Cro2;
            c(20) = CI2*Cro2;
            c(21) = CI2*Cro2;
            c(22) = CI2*Cro2;
            c(23) = CI2*RNAP;
            c(24) = CI2*RNAP;
            c(25) = CI2*RNAP;
            c(26) = Cro2*RNAP;
            c(27) = Cro2*RNAP;
            c(28) = Cro2*RNAP;
            c(29) = CI2^3;
            c(30) = Cro2^3;
            c(31) = CI2^2*Cro2;
            c(32) = CI2^2*Cro2;
            c(33) = CI2^2*Cro2;
            c(34) = CI2*Cro2^2;
            c(35) = CI2*Cro2^2;
            c(36) = CI2*Cro2^2;
            c(37) = RNAP*CI2^2;
            c(38) = RNAP*Cro2^2;
            c(39) = RNAP*CI2*Cro2;
            c(40) = RNAP*CI2*Cro2;
            
            weights = pg_P_RRM.*c;
            Z = sum(weights);
            probs = weights / Z;
        end
    end



end

function ind = pick_rxn(probs)
% Vectorized function to pick a state to jump to based on a list of
% probabilities

selection = rand;

M = cumsum([0, probs]);
M(end) = 1;

C = M >= selection;
D = M < selection;

ind = find(~(C(2:end) - D(1:end-1)));
end

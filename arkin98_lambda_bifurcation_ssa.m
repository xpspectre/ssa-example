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

% Avogadro's number
nA = 6.022e23; % (count/mol)

% Main states: x = [Cro, CI, CII, CIII, P1, P2, N, Cro2, CI2, P1_CII, P1_CIII, P2_CII, P2_CIII]
% Promotor states:
% Transcribing states:
% Translating states:

%% Initial conditions - counts
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

c0 = [Cro_0, CI_0, CII_0, CIII_0, P1_0, P2_0, N_0, Cro2_0, CI2_0, P1_CII_0, P1_CIII_0, P2_CII_0, P2_CIII_0];

RNAP_0 = 30e-9; % (M)
ribo_0 = 500e-9; % (M) ribosomes

x0 = round(c0*nA*V);

%% Promoter state
% P_RE: Table 1 of Arkin 98


% P_L: Table 1 of Arkin 98


% P_R and P_RM: Table 2 of Shea 85


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
            
            % Calculate volume
            %   Cell grows continuously from t = 0 to t = 35 min, doubling in size
            V = (1 + k0*t) * 1e-15; % (L)
            
            % Calculate promoter state probabilities
            
            
            % Calculate reaction propensities and normalization
            a = zeros(1,19); % total fixed reactions, additional ones based on txtl will be appended as needed
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

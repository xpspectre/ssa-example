function bistable_competition
% Positive feedback and mutual inhibition of 2 genes causes bistability
% Everything in normalized units
% Reactions:
%   0 -> A with rate k/(1+B^2)
%   A -> 0 with rate 1
%   0 -> B with rate k/(1+A^2)
%   B -> 0 with rate 1
% See T.S.Gardner et al, Nature 2000. From Problem Set 2, Problem 3b from
%   Systems Bio 7.32, 2016
clear; close all; clc
rng('default'); % remove this to get different traces

% Rate constants
%   k = 100 ensures that 1 will dominate forever (more or less); set tf >= 20
%   k = 10 allows switching between which one dominates; set tf >= 100
%   k = 1 means you can't see anything; set tf >= 1000
k = 100;

% Initial conditions - concs
% Start them both at the same amount and see if 1 dominates
%   This is completely fair but makes the ODE lines harder to distinguish
%   More interesting: set them slightly different and see the "wrong" one
%       dominate in the stochastic sim
A0 = 52;
B0 = 48;
x0 = [A0, B0];

% Set simulation conditions
tf = 30;

%% Run ODE simulation
[tOde, yOde] = ode15s(@ode_model, [0, tf], x0');

%% Run a few stochastic simulation
%   Run enough that (w/ a reproducible RNG seed) we see both A and B dominant runs
nSims = 3;
for iSim = 1:nSims
    [tSsa, xSsa] = ssa_sim(tf, x0);
    
    % Plot trajectories
    figure
    plot(tOde, yOde, 'LineWidth', 2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    stairs(tSsa, xSsa)
    hold off
    xlabel('Time')
    ylabel('Counts')
    title(sprintf('Bistable Competition, Slightly Off Initial Amounts\nODE = thick, SSA = thin lines'))
    xlim([0, tf])
    legend('A','B', 'Location','best')
end

%% Run sims with exactly matched initial A and B
% It's hard to see the A and B ODE lines exactly on top of each other
A0 = 50;
B0 = 50;
x0 = [A0, B0];

[tOde, yOde] = ode15s(@ode_model, [0, tf], x0');

nSims = 3;
for iSim = 1:nSims
    [tSsa, xSsa] = ssa_sim(tf, x0);
    
    % Plot trajectories
    figure
    plot(tOde, yOde, 'LineWidth', 2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    stairs(tSsa, xSsa)
    hold off
    xlabel('Time')
    ylabel('Counts')
    title(sprintf('Bistable Competition, Equal Initial Amounts\nODE = thick, SSA = thin lines'))
    xlim([0, tf])
    legend('A','B', 'Location','best')
end

% TODO: distributions of results

    function [ts, xs] = ssa_sim(tf, x0)
        % Gillespie stochastic simulation algorithm
        % States: x1 = A, x2 = B
        nt = 1e3; % preallocate times to start with
        ts = zeros(nt,1);
        it = 1;
        t = 0;
        ts(it) = t;
        
        nx = length(x0);
        xs = zeros(nt,nx); % preallocate solution
        x = x0;
        xs(it,:) = x;
        
        while t <= tf
            % Increment solution index
            it = it + 1;
            
            % Reallocate more solution space if needed (double each time)
            if it == nt
                ts = [ts; zeros(nt,1)];
                xs = [xs; zeros(nt,nx)];
                nt = nt * 2;
            end
            
            % Assign species to convenient vars
            A = x(1);
            B = x(2);
            
            % Calculate reaction propensities and normalization
            a = zeros(1,4);
            a(1) = k/(1 + B.^2);
            a(2) = A;
            a(3) = k/(1 + A.^2);
            a(4) = B;
            ao = sum(a);
            
            % Calculate next reaction time
            tau = (1/ao)*log(1/rand());
            
            % Update time
            t = t + tau;
            
            % Pick reaction
            u = pick_rxn(a/ao);
            
            % Handle corner case of no rxns possible (i.e., depletion)
            if isinf(t) || isempty(u)
                ts(it) = tf;
                xs(it,:) = x;
                break
            end
            
            % Update species
            switch u
                case 1 % 0 -> A
                    A = A + 1;
                case 2 % A -> 0
                    A = A - 1;
                case 3 % 0 -> B
                    B = B + 1;
                case 4 % B -> 0
                    B = B - 1;
            end
            
            % Reassign convenient vars to species
            x = [A, B];
            
            % Store results
            ts(it) = t;
            xs(it,:) = x;
        end
        
        % Trim solution
        ts = ts(1:it);
        xs = xs(1:it,:);
    end

    function dx = ode_model(~, x)
        % ODE model for integrators
        A = x(1);
        B = x(2);
        dA = k./(1 + B.^2) - A;
        dB = k./(1 + A.^2) - B;
        dx = [dA; dB];
    end

end

function state = pick_rxn(probs)
% Vectorized function to pick a state to jump to based on a list of
% probabilities

selection = rand;

M = cumsum([0, probs]);
M(end) = 1;

C = M >= selection;
D = M < selection;

state = find(~(C(2:end) - D(1:end-1)));
end

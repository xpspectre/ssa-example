function gardner00_bistable_competition
% Positive feedback and mutual inhibition of 2 genes causes bistability
%   based on the original syn bio paper
% Everything in normalized units, same forms as the paper, including
%   IPTG-tuned control
% Reactions:
%   0 -> u with rate alpha1/(1 + v^beta)
%   u -> 0 with rate 1
%   0 -> v with rate alpha2/(1 + u^gamma) or alpha2/(1 + (u/(1+IPTG/K)^eta)^gamma)
%   v -> 0 with rate 1
% See Gardner, T. S., Cantor, C. R., & Collins, J. J. (2000). Construction
%   of a genetic toggle switch in Escherichia coli. Nature, 403(6767), 339–342. http://doi.org/10.1038/35002131 
clear; close all; clc
rng('default'); % remove this to get different traces

% Rate constants
alpha1 = 156.25;
alpha2 = 15.6;
beta = 2.5;
gamma = 1;
eta = 2.0015;
K = 2.9618e-5;

IPTG = 0; % this is also basically a param since it's a constant in any expt

% Initial conditions - concs
u0 = 50;
v0 = 50;
x0 = [u0, v0];

% Set simulation conditions
tf = 24; % hr

%% Run ODE simulation
[tOde, yOde] = ode15s(@ode_model, [0, tf], x0');

% Get max val for plotting consistency
yMax = max(yOde(:));

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
    ylim([0, yMax*1.3])
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
        % States: x1 = u, x2 = v
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
            u = x(1);
            v = x(2);
            
            % Calculate reaction propensities and normalization
            a = zeros(1,4);
            a(1) = alpha1./(1 + v.^beta);
            a(2) = u;
            a(3) = alpha2./(1 + (u./(1+IPTG/K)^eta).^gamma);
            a(4) = v;
            ao = sum(a);
            
            % Calculate next reaction time
            tau = (1/ao)*log(1/rand());
            
            % Update time
            t = t + tau;
            
            % Pick reaction
            ind = pick_rxn(a/ao);
            
            % Handle corner case of no rxns possible (i.e., depletion)
            if isinf(t) || isempty(ind)
                ts(it) = tf;
                xs(it,:) = x;
                break
            end
            
            % Update species
            switch ind
                case 1 % 0 -> u
                    u = u + 1;
                case 2 % u -> 0
                    u = u - 1;
                case 3 % 0 -> v
                    v = v + 1;
                case 4 % v -> 0
                    v = v - 1;
            end
            
            % Reassign convenient vars to species
            x = [u, v];
            
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
        u = x(1);
        v = x(2);
        du = alpha1./(1 + v.^beta) - u;
        dv = alpha2./(1 + (u./(1+IPTG/K)^eta).^gamma) - v;
        dx = [du; dv];
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

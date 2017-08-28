function degradation_ssa
% Even simpler 1st order degradation rxn
% A -> 0
clear; close all; clc
rng('default');

% Rate constants
k = 5;

% Compartment volume - doesn't matter in this model
V = 1;

% Initial conditions - counts
%   Smaller for more stochasticity, larger for less
%   Timescale will always be the same
A0 = 50;

% Set simulation conditions
x0 = A0;
tf = 1;

% Run ODE simulation
[tOde, yOde] = ode15s(@ode_model, [0, tf], x0);

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
plot(tOde, yOde, 'LineWidth', 2)
hold on
ax = gca;
for i = 1:nRuns
    ax.ColorOrderIndex = 1;
    stairs(ts{i}, xs{i}) % using regular plot here is not as accurate
end
hold off
xlabel('Time')
ylabel('Counts')
title(sprintf('A -> 0\nODE = thick, SSA = thin lines'))
xlim([0, tf])
legend('A','Location','best')


    function [ts, xs] = ssa_sim(tf, x0)
        % Gillespie stochastic simulation algorithm
        % States: x1 = A
        % Rxns: a1 = degrade
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
            
            % Calculate reaction propensities and normalization
            a = k*x(1);
            ao = sum(a);
            
            % Calculate next reaction time
            tau = (1/ao)*log(1/rand());
            
            % Update time
            t = t + tau;
            
            % Pick reaction
            u = pick_rxn(a/ao);
            
            % Update species
            switch u
                case 1 % A -> 0
                    x(1) = x(1) - 1;
            end
            
            % Store results
            ts(it) = t;
            xs(:,it) = x;
        end
        
        % Trim solution
        ts = ts(1:it);
        xs = xs(:,1:it);
    end

    function dx = ode_model(~, x)
        % ODE model for integrators
        dx = -k*x;
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

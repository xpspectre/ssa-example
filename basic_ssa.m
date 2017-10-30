function basic_ssa
% Basic reversible binding model, i.e., protein + ligand <-> complex
% P + L <-> C
clear; close all; clc
rng('default');

% Rate constants
kon = 1e5; % (1/(M*s))
koff = 1e-2; % (1/s)

% Compartment volume
V = 1e-15; % (L) ~E. coli cell volume

% Avogadro's number
nA = 6.022e23; % count/mol

% Initial conditions - counts
P0 = 3e-7; % (M)
L0 = 1e-6; % (M)
C0 = 0; % (M)

c0 = [P0, L0, C0]; % in concs
x0 = round(c0*nA*V); % in counts

% Set simulation conditions
tf = 100;

% Run ODE simulation
[tOde, cOde] = ode15s(@ode_model, [0, tf], c0');
yOde = cOde*nA*V;

% Run stochastic simulation
%   Each run implicitly uses a different RNG seed
nRuns = 5;
ts = cell(nRuns,1);
xs = cell(nRuns,1);
for i = 1:nRuns
    [tSsa, xSsa] = ssa_sim(tf, x0);
    ts{i} = tSsa;
    xs{i} = xSsa;
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
title(sprintf('P + L <-> C\nODE = thick, SSA = thin lines'))
xlim([0, tf])
legend('P','L','C','Location','best')


    function [ts, xs] = ssa_sim(tf, x0)
        % Gillespie stochastic simulation algorithm
        % States: x1 = P, x2 = L, x3 = C
        % Rxns: a1 = on, a2 = off
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
%                 fprintf('Extended solution in step %i\n', it)
            end
            
            % Calculate reaction propensities and normalization
            a = zeros(1,2);
            a(1) = kon/(nA*V)*x(1)*x(2);
            a(2) = koff*x(3);
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
                xs(it,:) = x;
                break
            end
            
            % Update species
            switch u
                case 1 % A + B -> C
                    x(1) = x(1) - 1;
                    x(2) = x(2) - 1;
                    x(3) = x(3) + 1;
                case 2 % C -> A + B
                    x(1) = x(1) + 1;
                    x(2) = x(2) + 1;
                    x(3) = x(3) - 1;
            end
            
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
        dx1 = -kon*x(1)*x(2) + koff*x(3);
        dx2 = -kon*x(1)*x(2) + koff*x(3);
        dx3 =  kon*x(1)*x(2) - koff*x(3);
        dx = [dx1; dx2; dx3];
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

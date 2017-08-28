function basic_ssa
% Basic reversible binding model
% A + B <-> C
clear; close all; clc
rng('default');

% Rate constants
kon = 5;
koff = 3;

% Compartment volume
V = 1.5;

% Initial conditions - counts
A0 = 50;
B0 = 25;
C0 = 0;

% Set simulation conditions
x0 = [A0; B0; C0];
tf = 5e-2;

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
title(sprintf('A + B <-> C\nODE = thick, SSA = thin lines'))
xlim([0, tf])
legend('A','B','C','Location','best')


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
%                 fprintf('Extended solution in step %i\n', it)
            end
            
            % Calculate reaction propensities and normalization
            a = zeros(1,2);
            a(1) = kon*x(1)*x(2)/V;
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
                xs(:,it) = x;
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
            xs(:,it) = x;
        end
        
        % Trim solution
        ts = ts(1:it);
        xs = xs(:,1:it);
    end

    function dx = ode_model(~, x)
        % ODE model for integrators
        dx1 = -kon/V*x(1)*x(2) + koff*x(3);
        dx2 = -kon/V*x(1)*x(2) + koff*x(3);
        dx3 =  kon/V*x(1)*x(2) - koff*x(3);
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

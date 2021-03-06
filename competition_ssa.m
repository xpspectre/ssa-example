function competition_ssa
% Competition between 2 irreversible association reactions
% A + A -> AA
% A + B -> AB
% Handles some corner cases in the ODE eqs and SSA implementation
clear; close all; clc
rng('default');

% Rate constants
kAA = 1e6; % (1/(M*s))
kAB = 2e6; % (1/(M*s))

% Compartment volume
V = 1e-15; % (L) ~E. coli cell volume

% Avogadro's number
nA = 6.022e23; % count/mol

% Initial conditions - concs
A0 = 2e-7; % (M = mol/L)
B0 = 1e-7; % (M)
AA0 = 0;
AB0 = 0;

c0 = [A0; B0; AA0; AB0]; % in concs
x0 = round(c0*nA*V); % in counts
% 1e-6 M at 1e-15 L is ~ 600 molecules

% Set simulation conditions
tf = 60; % (s)

% Run ODE simulation
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
[tOde, cOde] = ode15s(@ode_model, [0, tf], c0, opts); % in concs
yOde = cOde*nA*V; % in counts

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
title(sprintf('A + A -> AA, A + B -> AB\nODE = thick, SSA = thin lines'))
xlim([0, tf])
legend('A','B','AA','AB', 'Location','best')

% Ensure mass balance is satisifed at all times
AOde = yOde(:,1) + 2*yOde(:,3) + yOde(:,4);
BOde = yOde(:,2) + yOde(:,4);

figure
plot(tOde, [AOde, BOde], 'LineWidth', 2)
hold on
ax = gca;
for i = 1:nRuns
    xSsa = xs{i};
    ASsa = xSsa(:,1) + 2*xSsa(:,3) + xSsa(:,4);
    BSsa = xSsa(:,2) + xSsa(:,4);
    ax.ColorOrderIndex = 1;
    stairs(ts{i}, [ASsa, BSsa], 'd') % using regular plot here is not as accurate
end
hold off
xlabel('Time')
ylabel('Counts')
title(sprintf('A + A -> AA, A + B -> AB\nMass Balance Validation'))
xlim([0, tf])
legend('A_{tot} ODE','B_{tot} ODE','A_{tot} SSA','B_{tot} SSA', 'Location','best')


    function [ts, xs] = ssa_sim(tf, x0)
        % Gillespie stochastic simulation algorithm
        % States: x1 = A, x2 = B, x3 = AA, x4 = AB
        % Rxns: a1 = AA rxn, a2 = AB rxn
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
            a = zeros(1,2);
            a(1) = kAA/(nA*V)*x(1)*(x(1)-1); % the 2*AA and x(1)*(x(1)-1)/2 cancels out
            a(2) = kAB/(nA*V)*x(1)*x(2);
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
                xs(:,it) = x;
                break
            end
            
            % Update species
            switch u
                case 1 % A + A -> AA
                    x(1) = x(1) - 2;
                    x(3) = x(3) + 1;
                case 2 % A + B -> AB
                    x(1) = x(1) - 1;
                    x(2) = x(2) - 1;
                    x(4) = x(4) + 1;
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
        dx1 = -2*kAA*x(1).^2 - kAB*x(1).*x(2); % A
        dx2 = -kAB*x(1).*x(2); % B
        dx3 =  kAA*x(1).^2; % AA
        dx4 =  kAB*x(1).*x(2); % AB
        dx = [dx1; dx2; dx3; dx4];
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

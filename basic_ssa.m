function basic_ssa
%% Basic reversible binding model, i.e., protein + ligand <-> complex
% P + L <-> C
clear; close all; clc
% rng('default');

% Rate constants
kon = 1e5; % (1/(M*s))
koff = 1e-2; % (1/s)

% Avogadro's number
nA = 6.022e23; % count/mol

% Initial conditions - counts
P0 = 3e-7; % (M)
L0 = 1e-6; % (M)
C0 = 0; % (M)

% Set simulation conditions
tf = 100; % (s)

% Compartment volume
% V = 1e-15; % (L) 1e-15 ~E. coli cell volume
% V = 1e-17;

%% Sims at different system sizes: same concs, different volumes
% Using same conc, different volume is easy (single param) and makes y-axes in
%   concs automatically the same
% OTOH, using same volume, different conc may be easier to conceptualize but
%   requires scaling y-axis and modifying all x0's.
useConcs = true; % if false, use counts

nRuns = 1000; % per volume
nShowTraces = 3;

% Vs = [1e-15, 1e-16, 1e-17];
% Vs = [1e-15, 5e-16, 1e-16];
Vs = [5e-17, 25e-17, 125e-17];
nV = length(Vs);

tOdes = cell(nV,1);
xOdes = cell(nV,1);
tSsass = cell(nV,1);
xSsass = cell(nV,1);
for iV = 1:nV
    V = Vs(iV);
    [tOde, xOde, tSsas, xSsas] = run_sims(V, nRuns);
    
    if useConcs
        xOde = xOde./(nA*V);
        for i = 1:nRuns
            xSsas{i} = xSsas{i}./(nA*V);
        end
    end
    
    tOdes{iV} = tOde;
    xOdes{iV} = xOde;
    tSsass{iV} = tSsas;
    xSsass{iV} = xSsas;
    
    figure
    plot(tOde, xOde, 'LineWidth', 2)
    hold on
    ax = gca;
    for i = 1:nShowTraces
        ax.ColorOrderIndex = 1;
        stairs(tSsas{i}, xSsas{i})
    end
    hold off
    xlabel('Time (s)')
    if useConcs
        ylabel('Conc (M)')
    else
        ylabel('Count')
    end
    title(sprintf('P + L <-> C, V = %.2g L\nCounts: P0 = %i, L0 = %i', V, round(P0*nA*V), round(L0*nA*V)))
    xlim([0, tf])
    legend('P','L','C','Location','best')
end

%% Overlaid +/- 1 std dev for each system size
% Interpolate at query points
nt = 100;
ts = linspace(0, tf, nt)';

figure
plot(tOde, xOde, 'LineWidth', 1)
hold on
ax = gca;
ymax = 0; % for setting the ylims
for iV = 1:nV
    xs = zeros(nRuns,nt,3);
    for i = 1:nRuns
        t = tSsass{iV}{i};
        x = xSsass{iV}{i};
        xs(i,:,:) = interp1(t, x, ts, 'previous');
    end
    
    xmean = squeeze(mean(xs, 1));
    xstd = squeeze(std(xs, 0, 1));
    nstds = 2;
    xlo = xmean - nstds*xstd;
    xhi = xmean + nstds*xstd;
    
    ymax = max(ymax, max(xhi(:))); % highest vals in the sim
    
    for i = 1:3
        color = ax.ColorOrder(i,:);
        fill([ts; flipud(ts)], [xlo(:,i); flipud(xhi(:,i))], color, 'EdgeColor', 'none');
        alpha(0.25) % overlaid to give differential shading
    end
    
end
hold off
xlabel('Time (s)')
ylabel('Conc (M)')
ylim([0, ymax*1.1])
title(sprintf('P + L <-> C\nODE and SSA with +/-%i std dev for different volumes', nstds))
xlim([0, tf])
legend('P','L','C','Location','best')

%% Helper functions
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

    function [tOde, xOde, tSsas, xSsas] = run_sims(V, nRuns)
        % Run deterministic and stochastic simulations with different system
        % sizes (but same concs)
        c0 = [P0, L0, C0]; % in concs
        x0 = round(c0*nA*V); % in counts
        
        % Run ODE simulation
        opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
        [tOde, cOde] = ode15s(@ode_model, [0, tf], c0', opts);
        xOde = cOde*nA*V;
        
        % Run stochastic simulation
        %   Each run implicitly uses a different RNG seed
        tSsas = cell(nRuns,1);
        xSsas = cell(nRuns,1);
        for iRun = 1:nRuns
            [tSsa, xSsa] = ssa_sim(tf, x0);
            tSsas{iRun} = tSsa;
            xSsas{iRun} = xSsa;
        end
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

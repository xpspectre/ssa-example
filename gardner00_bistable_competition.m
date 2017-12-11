function gardner00_bistable_competition
% Positive feedback and mutual inhibition of 2 genes causes bistability
%   based on the original syn bio paper
% Everything in normalized/nondimensionalized units, same forms as the paper, including IPTG-tuned control
% The system here has bistability and is one of the pTAK plasmids (nominally pTAK117 when using the Fig 5 parameter values)
%
% Figure 1 details:
%   Promoter 1: PLs1con
%   Repressor 1: cIts lambda repressor -> "v"
%   Inducer 1: Heat (denatures repressor 1)
%   Promoter 2: Ptrc-2
%   Repressor 2: LacI -> "u"
%   Inducer 2: IPTG (inhibits repressor 2)
%   Reporter: GFP, measure of amount of Repressor 1 and Promoter 2 activity
%
% Reactions:
%   0 -> u with rate alpha1/(1 + v^beta)
%   u -> 0 with rate 1
%   0 -> v with rate alpha2/(1 + u^gamma) or alpha2/(1 + (u/(1+IPTG/K)^eta)^gamma)
%   v -> 0 with rate 1
%
% See Gardner, T. S., Cantor, C. R., & Collins, J. J. (2000). Construction
%   of a genetic toggle switch in Escherichia coli. Nature, 403(6767), 339�342. http://doi.org/10.1038/35002131 
clear; close all; clc
rng('default'); % remove this to get different traces

% Rate constants
alpha1 = 156.25;
alpha2 = 15.6;
beta = 2.5;
gamma = 1;
eta = 2.0015;
K = 2.9618e-5; % (M = molar)

IPTG = 0; % (M) this is also basically a param since it's a constant in any expt

% Initial conditions
% u0 = 0;
% v0 = 0;
u0 = 50; % start off with something for example
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
    xlabel('Time (hr)')
    ylabel('Counts')
    title(sprintf('Bistable Competition\nODE = thick, SSA = thin lines'))
    xlim([0, tf])
    ylim([0, yMax*1.3])
    legend('u', 'v', 'Location','best')
end

%% Run sims with varying induction using IPTG
IPTGs = [0, 1e-6, 1e-2];
nIPTGs = length(IPTGs);

[tOde, yOde] = ode15s(@ode_model, [0, tf], x0');

for iIPTG = 1:nIPTGs
    IPTG = IPTGs(iIPTG);
    [tSsa, xSsa] = ssa_sim(tf, x0);
    
    % Plot trajectories
    figure
    plot(tOde, yOde, 'LineWidth', 2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    stairs(tSsa, xSsa)
    hold off
    xlabel('Time (br)')
    ylabel('Counts')
    title(sprintf('Bistable Competition, IPTG = %g\nODE = thick, SSA = thin lines', IPTG))
    xlim([0, tf])
    legend('u', 'v', 'Location','best')
end

%% Make Fig 5a
% Except instead of data, do multiple runs of the stochastic sim
% The lines are ODEs
u0 = 0; % there's nothing initially
v0 = 0;
x0 = [u0, v0];
tf = 17; % (hr) doesn't really make a difference compared to 24 hr

% Do a bunch of ODE sims to make a smooth curve (except at the jump)
nIPTGs = 50;
IPTGs = logspace(-6, -2, nIPTGs);
ssOdeOuts = nan(nIPTGs,1); % the reporter is just repressor 1, v
for iIPTG = 1:nIPTGs
    IPTG = IPTGs(iIPTG);
    [~, yOde] = ode15s(@ode_model, [0, tf], x0');
    ssOdeOuts(iIPTG) = yOde(end,2);
end

% Find jump for split plots
diffs = ssOdeOuts(2:end) - ssOdeOuts(1:end-1);
jump = find(diffs==max(diffs));

IPTGsLeft = IPTGs(1:jump);
IPTGsRight = IPTGs(jump+1:end);
ssOdeLeft = ssOdeOuts(1:jump);
ssOdeRight = ssOdeOuts(jump+1:end);

% Do a bunch of SSA sims at a few points and get stats
runSsas = false; % cache since this is expensive
ssasFile = 'gardner00_ssas_runs_.mat';

nIPTGsSsa = 10;
IPTGsSsa = logspace(-6, -2, nIPTGsSsa);
nSims = 100;

if runSsas
    ssSsaOuts = nan(nIPTGsSsa,nSims);
    for iIPTG = 1:nIPTGsSsa
        IPTG = IPTGsSsa(iIPTG);
        for iSim = 1:nSims
            [~, xSsa] = ssa_sim(tf, x0);
            ssSsaOuts(iIPTG,iSim) = xSsa(end,2);
        end
    end
    save(ssasFile, 'ssSsaOuts');
else
    loaded = load(ssasFile);
    ssSsaOuts = loaded.ssSsaOuts;
end

labels = cell(nIPTGsSsa,1);
for iIPTG = 1:nIPTGsSsa
    labels{iIPTG} = sprintf('%.2e', IPTGsSsa(iIPTG));
end

figure
boxplot(ssSsaOuts', 'Positions', log10(IPTGsSsa), 'PlotStyle', 'compact', 'Labels', labels) % median, IQR, limits, and outliers
hold on
ax = gca;
ax.ColorOrderIndex = 2; % red
plot(log10(IPTGsLeft), ssOdeLeft, 'LineWidth', 3)
ax.ColorOrderIndex = 2;
plot(log10(IPTGsRight), ssOdeRight, 'LineWidth', 3)
hold off
xlabel('[IPTG] (M)')
ylabel('cIts (counts)')
title(sprintf('Steady State Behavior\nODE = red line; SSA = blue boxplots'))

%% Distributions of steady state behaviors like Fig 5c
% Get domain for KDE
xLims = [0, max(ssSsaOuts(:))];
xSpan = xLims(2) - xLims(1);
xFactor = 0.25;
xLims = [xLims(1) - xFactor*xSpan, xLims(2) + xFactor*xSpan];

nxs = 100;
xs = linspace(xLims(1), xLims(2), nxs);

% Take a subset of IPTG concs since 10 is too many for the line plots
keep = 1:2:nIPTGsSsa;
nIPTGsSsa = length(keep);
ssSsaOuts = ssSsaOuts(keep,:);
labels = labels(keep);

fs = nan(nIPTGsSsa,nxs);
for iIPTG = 1:nIPTGsSsa
%     [fs(iIPTG,:),~] = ksdensity(ssSsaOuts(iIPTG,:), xs);
    [fs(iIPTG,:),~] = ksdensity(ssSsaOuts(iIPTG,:), xs, 'Bandwidth', 2); % these aren't actually normally distributed
end

colors = genColors(nIPTGsSsa);

figure
hold on
for iIPTG = 1:nIPTGsSsa
    plot(xs, fs(iIPTG,:), 'Color', colors(iIPTG,:))
end
hold off
xlabel('cIts (counts)')
ylabel('Prob')
h = legend(labels);
title(h, '[IPTG] (M)')
title('Distributions of Steady States')

%% Sims for phase portrait
% Modify alpha1 and alpha2 and do sims on each grid point
%   Do ODE sims to get deterministic behavior for basic plot - color indicates bistability?
%   Do many SSA sims to get distribution - mark each point by a color and shade?
% Show 2 regions bistability vs monostability
%   Bi: when you start at lo and hi, and get lo and hi results
%   Mono: when you start at lo and hi, and get 1 end behavior
na1 = 30;
na2 = 31;
alpha1s = logspace(log10(alpha1/100), log10(alpha1*100), na1);
alpha2s = logspace(log10(alpha2/100), log10(alpha2*100), na2);
IPTG = 0; % no induction - everything determined by starting state

vfs = nan(na1,na2,2); % last dim is [lo,hi]
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
for ia1 = 1:na1
    alpha1 = alpha1s(ia1);
    for ia2 = 1:na2
        alpha2 = alpha2s(ia2);
        
        % Common u start
        u0 = 0;
        
        % Lo cIts (v) start
        v0 = 0;
        [tOdeLo, yOdeLo] = ode15s(@ode_model, [0, tf], [u0, v0]', opts);
        vfs(ia1,ia2,1) = yOdeLo(end,2);
        
        % Hi cIts (v) start
        v0 = 100;
        [tOdeHi, yOdeHi] = ode15s(@ode_model, [0, tf], [u0, v0]', opts);
        vfs(ia1,ia2,2) = yOdeHi(end,2);
        
        % Plot the traces for specified cases
        % definite bistability
        if ia1 == na1 && ia2 == na2 % max vals for both - bistable
            figure
            plot(tOdeLo, yOdeLo(:,2), 'LineWidth', 2)
            hold on
            plot(tOdeHi, yOdeHi(:,2), 'LineWidth', 2)
            hold off
            xlabel('Time')
            ylabel('cIts')
            legend('Lo cIts Start','Hi cIts Start', 'Location','best')
            title('Hi \alpha_1, Hi \alpha_2 Traces')
        end
        
        if ia1 == na1 && ia2 == 1 % a corner - monostable
            figure
            plot(tOdeLo, yOdeLo(:,2), 'LineWidth', 2)
            hold on
            plot(tOdeHi, yOdeHi(:,2), 'LineWidth', 2)
            hold off
            xlabel('Time')
            ylabel('cIts')
            legend('Lo cIts Start','Hi cIts Start', 'Location','best')
            title('Hi \alpha_1, Lo \alpha_2 Traces')
        end
        
        if ia1 == 1 && ia2 == na2 % a corner - monostable
            figure
            plot(tOdeLo, yOdeLo(:,2), 'LineWidth', 2)
            hold on
            plot(tOdeHi, yOdeHi(:,2), 'LineWidth', 2)
            hold off
            xlabel('Time')
            ylabel('cIts')
            legend('Lo cIts Start','Hi cIts Start', 'Location','best')
            title('Lo \alpha_1, Hi \alpha_2 Traces')
        end
        
        if ia1 == 1 && ia2 == 1 % min vals for both - bistable?
            figure
            plot(tOdeLo, yOdeLo(:,2), 'LineWidth', 2)
            hold on
            plot(tOdeHi, yOdeHi(:,2), 'LineWidth', 2)
            hold off
            xlabel('Time')
            ylabel('cIts')
            legend('Lo cIts Start','Hi cIts Start', 'Location','best')
            title('Lo \alpha_1, Lo \alpha_2 Traces')
        end
    end
end

% Difference is a measure of bistability: same val means monostable;
%   "significant" diff means bistable
diffs = vfs(:,:,2) - vfs(:,:,1);

% Hack in a lower positive threshold for logspace plotting
eps = 1e-6;
diffs(diffs<0) = eps;

% Zscale in logscale
diffs = log10(diffs);

%%
figure
surf(alpha2s, alpha1s, diffs)
view(0,90);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
% set(gca, 'ZScale', 'log')
xlabel('\alpha_2')
ylabel('\alpha_1')
xlim([min(alpha2s), max(alpha2s)])
ylim([min(alpha1s), max(alpha1s)])
shading interp
h = colorbar;
ylabel(h, 'log_{10} Final cIts Diff')
title('Diff in Final cIts Between Starting Hi and Lo')

%% Output Figures
% Including setting proper overall figure size
% hs = findobj('Type', 'figure');
% for i = 1:8
%     fig = hs(i);
%     fig.PaperPositionMode = 'auto';
%     fig_pos = fig.PaperPosition;
%     fig.PaperSize = [fig_pos(3) fig_pos(4)];
%     saveas(fig, sprintf('fig%i.eps',i), 'epsc');
% end
for i = 1:13
    save2pdf(sprintf('fig%i.pdf',i), i);
end

%% End of main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Begin helper funs
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

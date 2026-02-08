function results = solve_spherical_cap_snapthrough()
%SOLVE_SPHERICAL_CAP_SNAPTHROUGH
% Axisymmetric spherical-cap BVP in the non-dimensional variables
% provided in the manuscript (no shallow-angle replacement in equations).
%
% Unknowns:
%   y = [w; w'; w''; Nphi; Ntheta; q]
% with
%   q = -sin(phi) w''' - cos(phi) w'' ...
%       + (nu sin(phi) + cos(phi)cot(phi)) w' ...
%       + Nphi w' sin(phi)
%
% Governing equations implemented exactly as written (in phi coordinate):
%   (Nphi sin phi)' - Ntheta cos phi = 0
%   q' = gamma (Nphi + Ntheta) sin phi         [point-load model]
%   tan(phi)*(Ntheta' - nu Nphi') + Ntheta(tan^2+nu+1)
%      - Nphi(nu tan^2+nu+1)
%      + 12(1-nu^2)*(0.5*(w')^2 - gamma*(w' tan + w tan^2)) = 0
%
% Apex load handling:
%   Default is true concentrated load via the apex boundary condition
%   lim_{phi->0} q = -gamma*Fbar/(2*pi), equivalent to the manuscript limit.
%   Optional patch regularization is kept as a switch.

    % -------- User parameters --------
    p.nu      = 0.30;      % Poisson ratio
    p.gamma   = 120;       % R/h
    p.phi0    = 0.35;      % cap semi-angle [rad]
    p.epsPole = 1e-4;      % start integration at small positive pole angle

    % Load model: 'point' (default, exact BC) or 'patch' (regularized).
    p.loadModel = 'point';
    p.phiPatch  = 0.02;    % used only when loadModel='patch'
    p.patchNormalization = 'domain'; % 'domain' uses [epsPole,phiPatch], 'fullcap' uses [0,phiPatch]
    p.patchLoadSign = +1;           % +1 or -1 for pressure sign convention in q' equation

    % Dimensionless total apex loads Fbar = F*R/D.
    FbarList = linspace(0, 120, 13);

    % Control mode: 'load' (prescribe FbarList) or 'displacement'
    % (prescribe apex displacement list and solve for Fbar).
    p.controlMode = 'load';
    p.wApexList = linspace(0, 2.5, 13); % used only in displacement mode
    p.forcePerturb = 1.0;               % finite-difference dF for displacement control
    p.dispTol = 1e-6;                   % tolerance on wbar(eps)-target
    p.maxDispNewton = 12;               % max Newton/secant iterations per target

    % Continuation / diagnostics controls (minimal stabilization).
    p.maxLoadStep = 10;      % automatic sub-stepping in Fbar
    p.residualNodes = 500;   % nodes for a posteriori residual estimate
    p.residualTrim = 0.02;   % trim fraction near each boundary for residual reporting
    p.residualQuantile = 0.95; % robust residual statistic (less spike-sensitive than max)

    % Mesh and initial guess.
    p.Fbar = FbarList(1);
    phiMesh = linspace(p.epsPole, p.phi0, 120);
    solinit = bvpinit(phiMesh, @(phi) guess_profile(phi, p));

    results = struct('Fbar', cell(numel(FbarList),1), ...
                     'sol',  cell(numel(FbarList),1), ...
                     'wApex',cell(numel(FbarList),1));

    opts = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',20000);

    switch lower(p.controlMode)
        case 'load'
            prevF = FbarList(1);
            for i = 1:numel(FbarList)
                targetF = FbarList(i);

                nSub = max(1, ceil(abs(targetF - prevF)/p.maxLoadStep));
                subLoads = linspace(prevF, targetF, nSub + 1);

                for k = 2:numel(subLoads)
                    p.Fbar = subLoads(k);

                    sol = bvp4c(@(phi,y) ode_system(phi,y,p), ...
                                @(ya,yb) bc_system(ya,yb,p), ...
                                solinit, opts);

                    solinit = bvpinit(sol.x, @(phi) deval(sol, phi));

                    yPoleStep = deval(sol, p.epsPole);
                    [odeResQ, odeResMax, bcResInf] = estimate_residuals(sol, p);
                    fprintf(['Step Fbar=%9.3f | wbar(eps)=%+12.6e | ', ...
                             'ODE-relres-q95=%9.3e | ODE-relres-max=%9.3e | BC-res-inf=%9.3e\n'], ...
                             p.Fbar, yPoleStep(1), odeResQ, odeResMax, bcResInf);
                end

                if targetF == prevF
                    p.Fbar = targetF;
                    sol = bvp4c(@(phi,y) ode_system(phi,y,p), ...
                                @(ya,yb) bc_system(ya,yb,p), ...
                                solinit, opts);
                    solinit = bvpinit(sol.x, @(phi) deval(sol, phi));
                    yPoleStep = deval(sol, p.epsPole);
                    [odeResQ, odeResMax, bcResInf] = estimate_residuals(sol, p);
                    fprintf(['Step Fbar=%9.3f | wbar(eps)=%+12.6e | ', ...
                             'ODE-relres-q95=%9.3e | ODE-relres-max=%9.3e | BC-res-inf=%9.3e\n'], ...
                             p.Fbar, yPoleStep(1), odeResQ, odeResMax, bcResInf);
                end

                yPole = deval(sol, p.epsPole);
                results(i).Fbar  = targetF;
                results(i).sol   = sol;
                results(i).wApex = yPole(1);

                prevF = targetF;
            end

        case 'displacement'
            nTar = numel(p.wApexList);
            results = struct('Fbar', cell(nTar,1), 'sol', cell(nTar,1), 'wApex', cell(nTar,1));

            Fguess = 0;
            for i = 1:nTar
                wTarget = p.wApexList(i);
                for it = 1:p.maxDispNewton
                    p.Fbar = Fguess;
                    sol = bvp4c(@(phi,y) ode_system(phi,y,p), ...
                                @(ya,yb) bc_system(ya,yb,p), ...
                                solinit, opts);
                    wNow = deval(sol, p.epsPole); wNow = wNow(1);
                    r = wNow - wTarget;

                    [odeResQ, odeResMax, bcResInf] = estimate_residuals(sol, p);
                    fprintf(['DispCtrl w*= %+10.4e | Fbar=%+10.4e | r=%+9.2e | ', ...
                             'ODE-relres-q95=%9.3e | ODE-relres-max=%9.3e | BC-res-inf=%9.3e\n'], ...
                             wTarget, Fguess, r, odeResQ, odeResMax, bcResInf);

                    if abs(r) < p.dispTol
                        break;
                    end

                    dF = max(p.forcePerturb, 1e-6*(1+abs(Fguess)));
                    p.Fbar = Fguess + dF;
                    solFD = bvp4c(@(phi,y) ode_system(phi,y,p), ...
                                  @(ya,yb) bc_system(ya,yb,p), ...
                                  bvpinit(sol.x, @(phi) deval(sol,phi)), opts);
                    wFD = deval(solFD, p.epsPole); wFD = wFD(1);
                    dwdF = (wFD - wNow)/dF;
                    if abs(dwdF) < 1e-12
                        dwdF = sign(dwdF + eps)*1e-12;
                    end

                    Fguess = Fguess - r/dwdF;
                    solinit = bvpinit(sol.x, @(phi) deval(sol,phi));
                end

                p.Fbar = Fguess;
                solinit = bvpinit(sol.x, @(phi) deval(sol,phi));
                results(i).Fbar = Fguess;
                results(i).sol = sol;
                results(i).wApex = wNow;
            end

        otherwise
            error('Unknown p.controlMode = %s. Use ''load'' or ''displacement''.', p.controlMode);
    end

    % Plot load-deflection curve.
    figure('Color','w');
    plot([results.Fbar], [results.wApex], 'o-','LineWidth',1.3);
    if strcmpi(p.controlMode, 'displacement')
        xlabel('Solved \mathcal{F}=FR/D');
        ylabel('Prescribed/achieved \bar{w}(\epsilon)');
    else
        xlabel('\mathcal{F}=FR/D');
        ylabel('\bar{w}(\epsilon) \approx \bar{w}(0)');
    end
    title(sprintf('Spherical-cap apex load-deflection (%s load model)', p.loadModel));
    grid on;

    % Example profile at highest load.
    figure('Color','w');
    phi = linspace(p.epsPole, p.phi0, 400);
    y   = deval(results(end).sol, phi);
    plot(phi, y(1,:), 'LineWidth',1.3);
    xlabel('\phi'); ylabel('\bar{w}(\phi)');
    title(sprintf('Deflection profile at \mathcal{F}=%.2f', results(end).Fbar));
    grid on;
end

function y = guess_profile(phi, p)
% Small-deflection clamped-cap-like guess.
    s = (phi - p.epsPole)/(p.phi0 - p.epsPole);
    w   = -0.02*(1 - s).^2;
    wp  =  0.04*(1 - s)/(p.phi0 - p.epsPole);
    wpp = -0.04/(p.phi0 - p.epsPole)^2;

    Np = 0.0;
    Nt = 0.0;

    if strcmpi(p.loadModel, 'point')
        q = -p.gamma*p.Fbar/(2*pi);
    else
        q = 0.0;
    end

    y = [w; wp; wpp; Np; Nt; q];
    y = y(:);
end

function dydphi = ode_system(phi, y, p)
% First-order system for y = [w; w'; w''; Nphi; Ntheta; q].
    nu    = p.nu;
    gamma = p.gamma;

    w   = y(1);
    wp  = y(2);
    wpp = y(3);
    Np  = y(4);
    Nt  = y(5);
    q   = y(6);

    s = sin(phi);
    c = cos(phi);
    t = tan(phi);
    cotp = c/s;

    % From q definition (exact trigonometric form).
    A = nu*s + c*cotp;
    wppp = (-q - c*wpp + A*wp + Np*wp*s)/s;

    % Membrane equilibrium.
    Np_p = (Nt - Np)*cotp;

    % Compatibility equation solved for Nt'.
    G = Nt*(t^2 + nu + 1) ...
      - Np*(nu*t^2 + nu + 1) ...
      + 12*(1 - nu^2)*(0.5*wp^2 - gamma*(wp*t + w*t^2));
    Nt_p = nu*Np_p - G/t;

    % Transverse balance.
    switch lower(p.loadModel)
        case 'point'
            q_p = gamma*(Np + Nt)*s;
        case 'patch'
            if phi <= p.phiPatch
                switch lower(p.patchNormalization)
                    case 'domain'
                        % Match total patch load on the *solved* interval [epsPole, phiPatch].
                        patchAreaFactor = cos(p.epsPole) - cos(p.phiPatch);
                    case 'fullcap'
                        % Match total patch load on [0, phiPatch] (legacy behavior).
                        patchAreaFactor = 1 - cos(p.phiPatch);
                    otherwise
                        error('Unknown p.patchNormalization = %s.', p.patchNormalization);
                end
                p0 = p.Fbar/(2*pi*patchAreaFactor);
            else
                p0 = 0;
            end
            q_p = gamma*(Np + Nt)*s + p.patchLoadSign*p0*s;
        otherwise
            error('Unknown p.loadModel = %s. Use ''point'' or ''patch''.', p.loadModel);
    end

    dydphi = [wp; wpp; wppp; Np_p; Nt_p; q_p];
end

function res = bc_system(ya, yb, p)
% BC residuals for 6 first-order equations.
% ya: at phi=epsPole, yb: at phi=phi0
    nu = p.nu;

    switch lower(p.loadModel)
        case 'point'
            qPoleTarget = -p.gamma*p.Fbar/(2*pi);
        case 'patch'
            switch lower(p.patchNormalization)
                case 'domain'
                    % No omitted load below epsPole in this normalization.
                    qPoleTarget = 0;
                case 'fullcap'
                    % Account for omitted 0..epsPole core load as equivalent concentrated jump.
                    omittedCore = 2*pi*(1 - cos(p.epsPole));
                    fullCapArea = 2*pi*(1 - cos(p.phiPatch));
                    Fcore = p.Fbar * omittedCore/fullCapArea;
                    qPoleTarget = -p.gamma*Fcore/(2*pi);
                otherwise
                    error('Unknown p.patchNormalization = %s.', p.patchNormalization);
            end
        otherwise
            error('Unknown p.loadModel = %s. Use ''point'' or ''patch''.', p.loadModel);
    end

    res = [ ...
        ya(2);                    % w'(eps)=0
        ya(4) - ya(5);            % Nphi=Ntheta at pole
        ya(6) - qPoleTarget;      % apex load BC in q-form
        yb(1);                    % w(phi0)=0
        yb(2);                    % w'(phi0)=0
        yb(5) - nu*yb(4)          % Ntheta=nu*Nphi at clamped edge
    ];
end


function [odeResQ, odeResMax, bcResInf] = estimate_residuals(sol, p)
%ESTIMATE_RESIDUALS A posteriori consistency check on converged BVP solution.
% Reports a relative ODE residual on a trimmed interior region to avoid
% singular-coordinate amplification near phi=epsPole.
    phi = linspace(p.epsPole, p.phi0, p.residualNodes);
    y = deval(sol, phi);

    dydphi_fd = zeros(size(y));
    for j = 1:size(y,1)
        dydphi_fd(j,:) = gradient(y(j,:), phi);
    end

    f = zeros(size(y));
    for i = 1:numel(phi)
        f(:,i) = ode_system(phi(i), y(:,i), p);
    end

    % Relative residual (component-wise), trimmed near boundaries.
    rel = abs(dydphi_fd - f)./(1 + abs(f));
    i0 = 1 + floor(p.residualTrim*(numel(phi)-1));
    i1 = numel(phi) - floor(p.residualTrim*(numel(phi)-1));
    i0 = max(i0, 2);
    i1 = max(i1, i0);
    relv = rel(:, i0:i1);
    relv = relv(:);
    odeResMax = max(relv);
    relv = sort(relv);
    kq = max(1, min(numel(relv), round(p.residualQuantile*numel(relv))));
    odeResQ = relv(kq);

    ya = y(:,1);
    yb = y(:,end);
    bc = bc_system(ya, yb, p);
    bcResInf = max(abs(bc));
end

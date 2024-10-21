%% File Setup Parameters & State Variables

%--General Parameters--%
%par.r = 0.0417;
par.r = .015;
par.R = 1+par.r;                                                           %Return on safe asset
par.beta = 0.97;                                                           %Discount factor in agents' utility function

par.penreplace = 0.75;                                                       %Pension Replacement Rate
par.smooth      = 2;                                                       

%-- Households --%
par.gamma       =5;                                                        % coefficient  of relative risk aversion in agents' utility function
par.U           = @(c) (c.^(1-par.gamma)-1)./(1-par.gamma);                % utility function
par.Uprime      = @(c) c.^(-par.gamma);                                    % first derivative of utility function
par.Uprimeinv   = @(u) u.^(-1./par.gamma);                                 % inverse of first derivative of utility function


%--Demographics--%
par.Twork       = 10;
par.Told        = 5;
par.T           = par.Twork+par.Told;

par.survprob    = ones(par.T,1);
par.survprob(par.T) = 0;
par.survprob(1:par.Twork) = linspace(1, 0.95, par.Twork);
par.survprob(par.Twork+1:par.T-1) =linspace(0.95, 0.55, par.T - par.Twork-1);


%-- Mean income profile --%
par.a = 0.045 * (10 / par.Twork);
par.b = 0.002 * (10/par.Twork)^2;

par.kappa    = par.a *[1:par.Twork]-par.b *([1:par.Twork].^2);
par.w = 1;                                                                 %Baseline wage


%--Asset Grid--%
par.amax = 30;
par.amin = 0;
par.Na = 25;
par.agrid =linspace(0,1,par.Na)';
agrid_par       = 0.4; %1 for linear, 0 for L-shaped
agrid           = par.agrid.^(1./agrid_par);
par.agrid       = par.amin + (par.amax-par.amin).*agrid;


%--Risky Asset Grid--%
par.risk_amin = 0;                                                         
par.risk_amax = 30;                                                        
par.Na_risk = 25;                                                          
par.risk_agrid = linspace(0, 1, par.Na_risk)';                         
agrid_risk_par = 0.4; % 1 for linear, 0 for L-shaped
agrid_risk = par.risk_agrid.^(1./agrid_risk_par);
par.risk_agrid = par.risk_amin + (par.risk_amax - par.risk_amin) .* agrid_risk;


%--Matrix A and A_prime-- %
%A is a (Na * Na_risk) * 2 matrix that contains each combination of assets
[safe_grid, risky_grid] = meshgrid(par.agrid, par.risk_agrid);
A = [safe_grid(:), risky_grid(:)];
A_prime = [safe_grid(:), risky_grid(:)];

%--Productivity Grid --%
par.Nz    = 5;
par.mu    = 0;                                                             % mean
par.rho   = 0.8;                                                           % persistence
par.sigma = 0.25;                                                          % st deviation
  
[par.zgrid, par.Pi] = tauchen(par);  

%--Risky Return Grid --%
par.Nr_risky    = 4;  
par.mu_risky    = par.R +0.04;                                             % mean
par.rho_risky   = 0;                                                       % persistence
par.sigma_risky = 0.2;                                                     % st deviation

[par.rgrid_risky, par.Pi_risk] = tauchen_risk(par);  

%par.Omega_risk = sparse(kron(par.Pi_risk, speye(par.Na * par.Na_risk)));
%par.Omega_risk = kron(speye(par.Na * par.Na_risk), par.Pi_risk);

%--Wages--%
par.wgrid     = exp(par.kappa)*par.w.*exp(par.zgrid/par.smooth);

%--Pension--%
par.pengrid = par.penreplace * exp(par.kappa(par.Twork)) * par.w * exp(par.zgrid/par.smooth);  % Pension for each productivity state (1D vector)




%% Simulation
par.Nsim        = 5000;
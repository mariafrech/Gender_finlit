%% File Setup Parameters & State Variables

%--General Parameters--%
par.r = 0.02;
par.R = 1+par.r;                                                           %Return on safe asset
par.beta = 0.96;                                                           %Discount factor in agents' utility function

par.penreplace = 0.68212;                                                       %Pension Replacement Rate
par.smooth      = 2;                                                       


%-- Households --%

par.gamma       =7;                                                        % coefficient  of relative risk aversion in agents' utility function
par.U           = @(c) (c.^(1-par.gamma)-1)./(1-par.gamma);                % utility function
par.Uprime      = @(c) c.^(-par.gamma);                                    % first derivative of utility function
par.Uprimeinv   = @(u) u.^(-1./par.gamma);                                 % inverse of first derivative of utility function

%par.gamma = 10;  % coefficient of relative risk aversion
%par.psi = 0.5;   % elasticity of intertemporal substitution (from life_cycle code)
%par.rho = 1 / par.psi;  % inverse of elasticity of intertemporal substitution

%par.U = @(c) ((c.^(1-par.rho) - 1) / (1-par.rho)).^(1/(1-par.gamma));

% First derivative of utility function
%par.Uprime = @(c) (((c.^(1-par.rho) - 1) / (1-par.rho)).^((1/(1-par.gamma))-(1))) .* (c.^(-par.rho));

% Inverse of first derivative of utility function
%par.Uprimeinv = @(u) ((1-par.rho) * (((1-par.gamma)*u).^((1-par.rho)/(1-par.gamma))) + 1).^(1/(1-par.rho));



%--Demographics--%
par.Twork       = 46;
par.Told        = 35;
par.T           = par.Twork+par.Told;

%par.survprob    = ones(par.T,1);
%par.survprob(par.T) = 0;
%par.survprob(1:par.Twork) = linspace(1, 0.95, par.Twork);
%par.survprob(par.Twork+1:par.T-1) =linspace(0.95, 0.55, par.T - par.Twork-1);

par.survprob = [
    0.99845; 0.99839; 0.99833; 0.9983; 0.99827; 0.99826; 0.99824; 0.9982; 
    0.99813; 0.99804; 0.99795; 0.99785; 0.99776; 0.99766; 0.99755; 0.99743; 
    0.9973; 0.99718; 0.99707; 0.99696; 0.99685; 0.99672; 0.99656; 0.99635; 
    0.9961; 0.99579; 0.99543; 0.99504; 0.99463; 0.9942; 0.9937; 0.99311; 
    0.99245; 0.99172; 0.99091; 0.99005; 0.98911; 0.98803; 0.9868; 0.98545; 
    0.98409; 0.9827; 0.98123; 0.97961; 0.97786; 0.97603; 0.97414; 0.97207; 
    0.9697; 0.96699; 0.96393; 0.96055; 0.9569; 0.9531; 0.94921; 0.94508; 
    0.94057; 0.9357; 0.93031; 0.92424; 0.91717; 0.90922; 0.90089; 0.89282; 
    0.88503; 0.87622; 0.86576; 0.8544; 0.8423; 0.82942; 0.8154; 0.80002; 
    0.78404; 0.76842; 0.75382; 0.73996; 0.72464; 0.71057; 0.6961; 0.6809; 0;
];



%{
%-- Mean income profile --%
par.a = 0.045 * (10 / par.Twork);
par.b = 0.002 * (10/par.Twork)^2;

par.kappa    = par.a *[1:par.Twork]-par.b *([1:par.Twork].^2);
par.w = 1;                                                                 %Baseline wage
%}

% Coefficients for the wage profile
aa = -2.170042 + 2.700381;  % constant term
b1 = 0.16818;               % linear coefficient
b2 = -0.0323371/10;         % quadratic coefficient
b3 = 0.0019704/100;         % cubic coefficient

par.kappa = aa + b1 *[20:par.Twork+20-1] + b2*([20:par.Twork+20-1].^2) + b3*([20:par.Twork+20-1].^3);
par.w = 1;

%--Asset Grid--%
par.amax = 200;
par.amin = 0.25;
par.Na = 50;
par.agrid =linspace(0,1,par.Na)';
%agrid_par       = 0.4; %1 for linear, 0 for L-shaped
%agrid           = par.agrid.^(1./agrid_par);
%par.agrid       = par.amin + (par.amax-par.amin).*agrid;
par.agrid = exp(linspace(log(par.amin), log(par.amax), par.Na))';          % Exponential grid


%--Risky Asset Grid--%
par.risk_amin = 0.25;                                                         
par.risk_amax = 200;                                                        
par.Na_risk = 50;                                                          
par.risk_agrid = linspace(0, 1, par.Na_risk)';                         
%agrid_risk_par = 0.4; % 1 for linear, 0 for L-shaped
%agrid_risk = par.risk_agrid.^(1./agrid_risk_par);
%par.risk_agrid = par.risk_amin + (par.risk_amax - par.risk_amin) .* agrid_risk;
par.risk_agrid = exp(linspace(log(par.risk_amin), log(par.risk_amax),par.Na_risk))';

%--Matrix A and A_prime-- %
%A is a (Na * Na_risk) * 2 matrix that contains each combination of assets
[safe_grid, risky_grid] = meshgrid(par.agrid, par.risk_agrid);
A = [safe_grid(:), risky_grid(:)];
A_prime = [safe_grid(:), risky_grid(:)];

% Set up A_prime
%[safe_grid, risky_grid] = ndgrid(par.agrid, par.risk_agrid);
%A = [safe_grid(:), risky_grid(:)];
%A_prime = [safe_grid(:), risky_grid(:)];

%--Productivity Grid --%
par.Nz    = 5;
par.mu    = 0;                                                             % mean
par.rho   = 0.8;                                                           % persistence
par.sigma = 0.25;                                                          % st deviation
  
[par.zgrid, par.Pi] = tauchen(par);  

%--Risky Return Grid --%
par.Nr_risky    = 3;  
par.mu_risky    = par.R +0.06;                                             % mean
par.rho_risky   = 0;                                                       % persistence
par.sigma_risky = 0.2;                                                     % st deviation

[par.rgrid_risky, par.Pi_risk] = tauchen_risk(par);  

%par.Omega_risk = sparse(kron(par.Pi_risk, speye(par.Na * par.Na_risk)));
par.Omega_risk = kron(speye(par.Na * par.Na_risk), par.Pi_risk);

%--Wages--%
par.wgrid     = exp(par.kappa + (par.zgrid/par.smooth)); % * 0.0825;
%par.wgrid = [1.0000, 1.1672, 1.3543, 1.5629, 1.7938, 2.0473, 2.3232, 2.6211, 2.9403, 3.2800, 3.6392, 4.0167, 4.4111, 4.8209, 5.2444, 5.6797, 6.1245, 6.5768, 7.0339, 7.4935, 7.9528, 8.4092, 8.8601, 9.3029, 9.7350, 10.1538, 10.5568, 10.9418, 11.3064, 11.6486, 11.9664, 12.2582, 12.5224, 12.7575, 12.9625, 13.1363, 13.2780, 13.3874, 13.4638, 13.5071, 13.5173, 13.4945, 13.4390, 13.3517, 13.2332, 13.0848];
%par.wgrid = exp(par.kappa)* par.w.* exp(par.zgrid/par.smooth);
%--Pension--%

par.pengrid = par.penreplace * exp(par.kappa(par.Twork)) * par.w * exp(par.zgrid/par.smooth);%* 0.0825;  % Pension for each productivity state (1D vector)

%par.pengrid = par.penreplace * par.wgrid;



%% Simulation
par.Nsim        = 5000;
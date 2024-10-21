function [rgrid_risky,Pi_risk] = tauchen_risk(par,bandwidth)
% This code computes the grid and transition matrix associated to an AR(1)
% process discretized according to Tauchen (1986). It approximates a
% the process
%       z(t+1) = (1-rho)*mu+rho*z(t)+eps(t+1) 
% where Var(eps)=sigma^2.
%
% Input:    N           scalar, number of gridpoints for z
%           mu          scalar, unconditional mean of the z process
%           rho         scalar, autocorrelation parameter
%           sigma       scalar, stdev of epsilon
%           bandwidth   positive integer, multiple of sigma to define bandwidth (optional, 3 by default)
%
% Output:   zgrid   Nx1 vector, grid for z
%           Pi      NxN matrix, transition probability matrix
%
% ------------------------------------------------------------------------
%             Advanced Techniques in Macroeconomics I
%                           -------------
% Instructor: Edouard Schaal
% ------------------------------------------------------------------------


if nargin<5 % if bandwidth is not specified
    bandwidth = 3; % set bandwidth to 3 stdev's by default
end

sigmaz = par.sigma_risky / sqrt(1-par.rho_risky^2); % stdev of the z process

% Create grid

rgrid_risky = linspace( par.mu_risky-bandwidth*sigmaz, par.mu_risky+bandwidth*sigmaz,par.Nr_risky)';

% Create transition matrix

brackets  = [-inf;0.5*(rgrid_risky(2:end) + rgrid_risky(1:end-1));inf]; % define the brackets that define the bins

[z_mat,brackets_mat]=ndgrid(rgrid_risky,brackets); % define two Nx(N+1) grids of z values and bracket values

Pi_risk = normcdf(brackets_mat(:,2:end) - par.rho_risky*z_mat(:,2:end) - (1-par.rho_risky)*par.mu_risky,0,par.sigma_risky)...
    -normcdf(brackets_mat(:,1:end-1) - par.rho_risky*z_mat(:,1:end-1) - (1-par.rho_risky)*par.mu_risky,0,par.sigma_risky);




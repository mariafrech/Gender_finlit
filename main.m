%--------------------------%
%Two-Asset Life-Cycle Model%
%--------------------------%

% ------------------------------------------------------------------------
% PREAMBLE
clear; clc; close all;
tic;

%% ------------------------------------------------------------------------
% OPTIONS
Display          = 1;
DoSimulate       = 1;
MakePlots        = 1;

% ------------------------------------------------------------------------
% PARAMETERS and FUNCTIONS
addpath('functions')
addpath('utils')
%param; 
param_cocco


%% ------------------------------------------------------------------------
% INITIALIZE
% Cell array for each exogenous state(prod, risky_return and time)

V_ret  = cell(par.Nz, par.Nr_risky, par.Told);                             % Value functions for retirement periods   
V_w    = cell(par.Nz, par.Nr_risky, par.Twork);                            % Value functions for working periods

sav_ret  = cell(par.Nz, par.Nr_risky, par.Told);                           % Savings policies for retirement
inv_ret  = cell(par.Nz, par.Nr_risky, par.Told);                           % Investment policies for retirement
cons_ret = cell(par.Nz, par.Nr_risky, par.Told);                           % Consumption during retirement periods

sav_w  = cell(par.Nz, par.Nr_risky, par.Twork);                            % Savings policies for working periods
inv_w  = cell(par.Nz, par.Nr_risky, par.Twork);                            % Investment policies for working periods
cons_w = cell(par.Nz, par.Nr_risky, par.Twork);                            % Consumption during working periods


% ------------------------------------------------------------------------
% Final period (T = par.Told)

for z = 1:par.Nz
    for r = 1:par.Nr_risky
        % Calculate consumption for all asset combinations
        cons_today = par.R * A(:,1) + par.rgrid_risky(r) * A(:,2) + par.pengrid(z);
        
       
        % Utility calculation
        V_final = par.U(cons_today);
        V_final(cons_today <= 0) = -inf;
        
        % Store results
        V_ret{z, r, par.Told} = V_final;
        cons_ret{z, r, par.Told} = cons_today;
        sav_ret{z, r, par.Told} = zeros(size(A(:,1)));
        inv_ret{z, r, par.Told} = zeros(size(A(:,2)));
    end
end

% ------------------------------------------------------------------------
%% Solve backwards for retirement periods

for t = par.Told-1:-1:1
    disp(['Solving retirement at age: ' int2str(par.Twork + t)]);

    
    for z = 1:par.Nz
        for r = 1:par.Nr_risky
            % Initialize value function for this state
            V_ret{z,r,t} = zeros(size(A,1), 1);
            sav_ret{z,r,t} = zeros(size(A,1), 1);
            inv_ret{z,r,t} = zeros(size(A,1), 1);
            
            % Loop over all possible current asset combinations
            for i = 1:size(A,1)
                a_today = A(i,:);
                
                % Calculate budget constraint
                budget = par.R * a_today(1) + par.rgrid_risky(r) * a_today(2) + par.pengrid(z);
               

                % Find optimal future asset allocation
               [v_max, idx] = max_future_value_ret(budget, A_prime, z, r, t, par, V_ret);
             % [v_max, sav_opt, inv_opt] = max_future_value_ret_interp(budget, A_prime, z, r, t, par, V_ret);
                
 
                % Store results
                V_ret{z,r,t}(i) = v_max;
                sav_ret{z,r,t}(i) = A(idx,1);
                inv_ret{z,r,t}(i) = A(idx,2);
                % sav_ret{z,r,t}(i) = sav_opt;
                %inv_ret{z,r,t}(i) = inv_opt;

              
            end
            

            % Calculate corresponding consumption
            cons_ret{z,r,t} = par.R * A(:,1) + par.rgrid_risky(r) * A(:,2) + par.pengrid(z) - sav_ret{z,r,t} - inv_ret{z,r,t};
            
          

        end
    end
end
   
plot_policy_functions_ret(sav_ret, inv_ret, cons_ret,par,34)

% Call the function after solving the model
%check_corner_solutions(sav_ret, inv_ret, A, par);

% After calculating policy functions
%check_budget_constraint(par, cons_ret, sav_ret, inv_ret, 1e-6);


%% ------------------------------------------------------------------------
%% Final working period (t = par.Twork)

disp(['Solving at working age: ' int2str(par.Twork)]);

for z = 1:par.Nz
    for r = 1:par.Nr_risky
        % Initialize value function for this state
        V_w{z,r,par.Twork} = zeros(size(A,1), 1);
        sav_w{z,r,par.Twork} = zeros(size(A,1), 1);
        inv_w{z,r,par.Twork} = zeros(size(A,1), 1);
        
        % Loop over all possible current asset combinations
        for i = 1:size(A,1)
            a_today = A(i,:);
            
            % Calculate budget constraint
            budget = par.R * a_today(1) + par.rgrid_risky(r) * a_today(2) + par.wgrid(z,par.Twork);
            
            % Find optimal future asset allocation
            [v_max, idx] = max_future_value_last_work(budget, A, z, r, par.Twork, par, V_ret);
            
            % Store results
            V_w{z,r,par.Twork}(i) = v_max;
            sav_w{z,r,par.Twork}(i) = A(idx,1);
            inv_w{z,r,par.Twork}(i) = A(idx,2);
        end
        
        % Calculate corresponding consumption
        cons_w{z,r,par.Twork} = par.R * A(:,1) + par.rgrid_risky(r) * A(:,2) + par.wgrid(z,par.Twork) - sav_w{z,r,par.Twork} - inv_w{z,r,par.Twork};
    end
end

%% ------------------------------------------------------------------------
%% Other working periods

for t = par.Twork-1:-1:1
        disp(['Solving at working age: ' int2str(t)]);
    
for z = 1:par.Nz
    for r = 1:par.Nr_risky
        % Initialize value function for this state
        V_w{z,r,t} = zeros(size(A,1), 1);
        sav_w{z,r,t} = zeros(size(A,1), 1);
        inv_w{z,r,t} = zeros(size(A,1), 1);
        
        % Loop over all possible current asset combinations
        for i = 1:size(A,1)
            a_today = A(i,:);
            
            % Calculate budget constraint
            budget = par.R * a_today(1) + par.rgrid_risky(r) * a_today(2) + par.wgrid(z,t);
            
            % Find optimal future asset allocation
            [v_max, idx] = max_future_value_work(budget, A, z, r, t, par, V_w);
            
            % Store results
            V_w{z,r,t}(i) = v_max;
            sav_w{z,r,t}(i) = A(idx,1);
            inv_w{z,r,t}(i) = A(idx,2);
        end
        
        % Calculate corresponding consumption
        cons_w{z,r,t} = par.R * A(:,1) + par.rgrid_risky(r) * A(:,2) + par.wgrid(z,t) - sav_w{z,r,t} - inv_w{z,r,t};
    end
end

end

%% ------------------------------------------------------------------------
% SAVE RESULTS

save('policy_functions.mat', 'par', 'cons_w', 'sav_w', 'inv_w', 'cons_ret', 'sav_ret', 'inv_ret', 'V_w', 'V_ret');
disp('Results have been saved to policy_functions.mat');

elapsed_time = toc;
fprintf('Elapsed time: %.2f seconds\n', elapsed_time);

if DoSimulate == 1
    [csim, asim, arsim, wagesim] = run_simulation_and_plot(par, sav_w, inv_w, cons_w, sav_ret, inv_ret, cons_ret, Display, MakePlots);
end




function check_budget_constraint(par, cons_ret, sav_ret, inv_ret, tolerance)
    [Na, Na_risk, Nz] = size(sav_ret);
    violations = zeros(Na, Na_risk, Nz);
    max_violations = zeros(Na, Na_risk, Nz);

    for i = 1:Na
        for j = 1:Na_risk
            for k = 1:Nz
                cons = cons_ret{i,j,k};
                sav = sav_ret{i,j,k};
                inv = inv_ret{i,j,k};
                
                pension = par.pengrid(k);
                
                left_side = cons + sav + inv;
                right_side = pension + (1 + par.r) * par.agrid(i) + (1 + par.rgrid_risky(j)) * par.risk_agrid(j);
                right_side = repmat(right_side, size(left_side));
                
                diff = abs(left_side - right_side);
                violations(i,j,k) = sum(diff > tolerance);
                max_violations(i,j,k) = max(diff);
            end
        end
    end
    
    [max_val, idx] = max(max_violations(:));
    [i_max, j_max, k_max] = ind2sub(size(max_violations), idx);
    
    fprintf('Maximum violation: %f at i=%d, j=%d, k=%d\n', max_val, i_max, j_max, k_max);
    fprintf('Total violations: %d\n', sum(violations(:)));
    
    % Plot heatmap of violations
    figure;
    imagesc(sum(violations, 3));
    colorbar;
    title('Heatmap of Budget Constraint Violations');
    xlabel('Risky Asset Index');
    ylabel('Safe Asset Index');
end



function check_corner_solutions(sav_ret, inv_ret, A, par)
    [Na, Na_risk, Nz, Nr_risky] = size(sav_ret);
    T = length(sav_ret{1,1,1,1});  % Assuming time is the last dimension in each cell
    
    for t = 1:T
        all_safe = 0;
        all_risky = 0;
        no_savings = 0;
        total_decisions = 0;
        
        for z = 1:Nz
            for r = 1:Nr_risky
                for i = 1:Na
                    for j = 1:Na_risk
                        total_wealth = A(i,1) + A(i,2);
                        safe_share = sav_ret{i,j,z,r}(t) / total_wealth;
                        risky_share = inv_ret{i,j,z,r}(t) / total_wealth;
                        
                        if safe_share > 0.99
                            all_safe = all_safe + 1;
                        elseif risky_share > 0.99
                            all_risky = all_risky + 1;
                        elseif safe_share + risky_share < 0.01
                            no_savings = no_savings + 1;
                        end
                        
                        total_decisions = total_decisions + 1;
                    end
                end
            end
        end
        
        fprintf('Time period %d:\n', t);
        fprintf('All safe assets: %.2f%%\n', 100 * all_safe / total_decisions);
        fprintf('All risky assets: %.2f%%\n', 100 * all_risky / total_decisions);
        fprintf('No savings: %.2f%%\n', 100 * no_savings / total_decisions);
        fprintf('\n');
        fprintf('State combination %d (Safe: %.2f, Risky: %.2f, z: %d, r_prev: %.2f):\n', ...
        i, safe_asset_value, risky_asset_value, productivity_state, previous_return);

    end
end


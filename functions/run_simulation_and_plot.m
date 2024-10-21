function [csim, asim, arsim, wagesim] = run_simulation_and_plot(par, sav_w, inv_w, cons_w, sav_ret, inv_ret, cons_ret, Display, MakePlots)
    % Run simulation for the life-cycle model and generate plots
    %
    % Inputs:
    %   par - Structure containing model parameters
    %   sav_w, inv_w, cons_w - Policy functions for working periods
    %   sav_ret, inv_ret, cons_ret - Policy functions for retirement periods
    %   Display - Flag to display progress messages
    %   MakePlots - Flag to generate plots
    %
    % Outputs:
    %   csim - Simulated consumption
    %   asim - Simulated safe assets
    %   arsim - Simulated risky assets
    %   wagesim - Simulated wages

    % Set random seed for reproducibility
    rng(2021);

    % Generate random initial wealth
    zrand = rand(par.Nsim, par.Twork);
    %arand = par.amin+1 + par.amax/100 * rand(par.Nsim, 1);
    %ariskrand = par.risk_amin + par.risk_amax/100 * rand(par.Nsim, 1);
    arand = zeros(par.Nsim, 1);
    ariskrand = zeros(par.Nsim, 1);

    % Preallocate matrices
    zindsim = zeros(par.Nsim, par.Twork);
    rindsim = zeros(par.Nsim, par.T);
    csim = zeros(par.Nsim, par.T);
    zsim = zeros(par.Nsim, par.Twork);
    asim = zeros(par.Nsim, par.T + 1);
    arsim = zeros(par.Nsim, par.T + 1);
    wagesim = zeros(par.Nsim, par.Twork);
    rsim = zeros(par.Nsim, par.T);

    % Initialize first period asset distribution with non-zero wealth
    asim(:, 1) = arand;  % Initial safe assets
    arsim(:, 1) = ariskrand;  % Initial risky assets

    % Generate productivity for each period for each 'agent'
    for n = 1:par.Nsim
        mc = dtmc(par.Pi);
        numSteps = par.Twork - 1;
        zindsim(n, :) = simulate(mc, numSteps);
        zsim(n, :) = exp(par.zgrid(zindsim(n, :)))';
    end

    % Generate risky return indices for each period for each 'agent'
    for n = 1:par.Nsim
        mc_risk = dtmc(par.Pi_risk);
        numSteps = par.T - 1;
        rindsim(n, :) = simulate(mc_risk, numSteps);
        rsim(n, :) = par.rgrid_risky(rindsim(n, :));
    end

    % Main simulation loop
    for it = 1:par.T
        if Display >= 1
            disp(['Simulating at age: ' int2str(it)]);
        end

        for n = 1:par.Nsim
            % Find the nearest grid points for current asset levels
            [~, a_idx] = min(abs(par.agrid - asim(n, it)));
            [~, ar_idx] = min(abs(par.risk_agrid - arsim(n, it)));

            % Combine safe and risky asset indices
            combined_idx = sub2ind([length(par.agrid), length(par.risk_agrid)], a_idx, ar_idx);

            if it <= par.Twork
                % Working period
                z_idx = zindsim(n, it);
                r_idx = rindsim(n, it);

                % Get optimal savings and investment
                asim(n, it+1) = sav_w{z_idx, r_idx, it}(combined_idx);
                arsim(n, it+1) = inv_w{z_idx, r_idx, it}(combined_idx);

                % Get consumption and wage
                csim(n, it) = cons_w{z_idx, r_idx, it}(combined_idx);
                wagesim(n, it) = par.wgrid(z_idx, it);
            else
                % Retirement period
                z_idx = zindsim(n, par.Twork);  % Use last working period productivity
                r_idx = rindsim(n, it);
                t_ret = it - par.Twork;  % Adjust index for retirement period

                % Get optimal savings and investment
                asim(n, it+1) = sav_ret{z_idx, r_idx, t_ret}(combined_idx);
                arsim(n, it+1) = inv_ret{z_idx, r_idx, t_ret}(combined_idx);

                % Get consumption
                csim(n, it) = cons_ret{z_idx, r_idx, t_ret}(combined_idx);
            end
        end
    end

    % Make plots
    if MakePlots == 1
        generate_plots(par, csim, asim, arsim, wagesim);
        generate_combined_plot(par, csim, asim, arsim, wagesim);
   
    end
end

function generate_plots(par, csim, asim, arsim, wagesim)
    figure(1);
    x = 1:par.T;

    % Consumption
    subplot(2, 2, 1);
    plot(x, mean(csim, 1), 'b--', 'LineWidth', 2);
    title('Consumption');
    xlabel('Age');
    ylabel('Mean Consumption');
    grid on;

    % Safe Assets
    subplot(2, 2, 2);
    plot(x, mean(asim(:, 2:end), 1), 'b--', 'LineWidth', 2);
    title('Safe Assets');
    xlabel('Age');
    ylabel('Mean Safe Assets');
    grid on;

    % Risky Assets
    subplot(2, 2, 3);
    plot(x, mean(arsim(:, 2:end), 1), 'b--', 'LineWidth', 2);
    title('Risky Assets');
    xlabel('Age');
    ylabel('Mean Risky Assets');
    grid on;

    % Wage
    subplot(2, 2, 4);
    plot(1:par.Twork, mean(wagesim, 1), 'r--', 'LineWidth', 2);
    title('Mean Wage');
    xlabel('Age');
    ylabel('Mean Wage');
    grid on;

    % Save figure
    output_folder = 'Figures';
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    saveas(gcf, fullfile(output_folder, 'simulation_results'));
    close(gcf);
end

function generate_combined_plot(par, csim, asim, arsim, wagesim)
    figure(2);
    x = 1:par.T;

    % Calculate total wealth (safe assets + risky assets)
    total_wealth = asim(:, 2:end) + arsim(:, 2:end);

    % Calculate wage including retirement (using mean wage of last working period * replacement rate)
    extended_wage = zeros(size(wagesim, 1), par.T);
    extended_wage(:, 1:par.Twork) = wagesim;
    mean_last_wage = mean(wagesim(:, par.Twork));
    extended_wage(:, par.Twork+1:end) = mean_last_wage * par.penreplace;

    % Plot combined graph
    plot(x, mean(csim, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'Consumption');
    hold on;
    plot(x, mean(total_wealth, 1), 'r-', 'LineWidth', 2, 'DisplayName', 'Total Wealth');
    plot(x, mean(extended_wage, 1), 'g-', 'LineWidth', 2, 'DisplayName', 'Wage/Pension');
    hold off;

    title('Combined Life-Cycle Plot');
    xlabel('Age');
    ylabel('Value');
    legend('show');
    grid on;

    % Save figure
    output_folder = 'Figures';
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    saveas(gcf, fullfile(output_folder, 'combined_life_cycle_plot'));
    close(gcf);
end


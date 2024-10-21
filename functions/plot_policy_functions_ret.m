%{
function plot_policy_functions_ret(sav_ret, inv_ret, cons_ret, par, t)
    % t is the time index to plot (e.g., par.Told - 10 for 10 periods before the end)
    
    % Find the median safe and risky asset indices
    median_safe_idx = round(par.Na / 2);
    median_risk_idx = round(par.Na_risk / 2);
    
    % Define colormap
    colors = jet(par.Nz * par.Nr_risky);
    
    % Create two figures
    figure_titles = {'Fixed Median Safe Asset', 'Fixed Median Risky Asset'};
    for fig = 1:2
        figure('Position', [100, 100, 1500, 400]);
        
        % Subplot for Savings
        subplot(1, 3, 1);
        hold on;
        title('Savings Policy');
        if fig == 1
            xlabel('Risky Asset');
        else
            xlabel('Safe Asset');
        end
        ylabel('Savings');
        
        % Subplot for Investment
        subplot(1, 3, 2);
        hold on;
        title('Investment Policy');
        if fig == 1
            xlabel('Risky Asset');
        else
            xlabel('Safe Asset');
        end
        ylabel('Investment');
        
        % Subplot for Consumption
        subplot(1, 3, 3);
        hold on;
        title('Consumption Policy');
        if fig == 1
            xlabel('Risky Asset');
        else
            xlabel('Safe Asset');
        end
        ylabel('Consumption');
        
        legend_labels = cell(par.Nz * par.Nr_risky, 1);
        color_idx = 1;
        
        for z = 1:par.Nz
            for r = 1:par.Nr_risky
                % Extract policy functions for this state
                sav = reshape(sav_ret{z, r, t}, par.Na, par.Na_risk);
                inv = reshape(inv_ret{z, r, t}, par.Na, par.Na_risk);
                cons = reshape(cons_ret{z, r, t}, par.Na, par.Na_risk);
                
                if fig == 1  % Fixed median safe asset
                    x_values = par.risk_agrid;
                    sav_plot = sav(median_safe_idx, :);
                    inv_plot = inv(median_safe_idx, :);
                    cons_plot = cons(median_safe_idx, :);
                else  % Fixed median risky asset
                    x_values = par.agrid;
                    sav_plot = sav(:, median_risk_idx);
                    inv_plot = inv(:, median_risk_idx);
                    cons_plot = cons(:, median_risk_idx);
                end
                
                % Plot Savings
                subplot(1, 3, 1);
                plot(x_values, sav_plot, 'Color', colors(color_idx, :), 'LineWidth', 2);
                
                % Plot Investment
                subplot(1, 3, 2);
                plot(x_values, inv_plot, 'Color', colors(color_idx, :), 'LineWidth', 2);
                
                % Plot Consumption
                subplot(1, 3, 3);
                plot(x_values, cons_plot, 'Color', colors(color_idx, :), 'LineWidth', 2);
                
                % Create legend label
                legend_labels{color_idx} = sprintf('z=%d, r=%d', z, r);
                color_idx = color_idx + 1;
            end
        end
        
        % Add legend to the last subplot
        subplot(1, 3, 3);
        legend(legend_labels, 'Location', 'eastoutside');
        
        % Adjust subplot layouts to make room for the legend
        set(gcf, 'Position', [100, 100, 1500, 400]);
        
        % Add an overall title
        sgtitle(sprintf('Policy Functions at t = %d, %s', par.Twork + t, figure_titles{fig}));
        
        % Save figure
        output_folder = 'Figures';
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        saveas(gcf, fullfile(output_folder, sprintf('policy_functions_%s', strrep(figure_titles{fig}, ' ', '_'))));
        close(gcf);
    end
end
%}

function plot_policy_functions_ret(sav_ret, inv_ret, cons_ret, par, t)
    % t is the time index to plot (e.g., par.Told - 10 for 10 periods before the end)
    
    % Find the median safe and risky asset indices
    median_safe_idx = round(par.Na / 2);
    median_risk_idx = round(par.Na_risk / 2);
    
    % Define colormap
    colors = jet(par.Nz);
    
    % Create two figures
    figure_titles = {'Fixed Median Safe Asset', 'Fixed Median Risky Asset'};
    for fig = 1:2
        figure('Position', [100, 100, 1500, 400]);
        
        % Subplot for Savings
        subplot(1, 3, 1);
        hold on;
        title('Savings Policy');
        if fig == 1
            xlabel('Risky Asset');
        else
            xlabel('Safe Asset');
        end
        ylabel('Savings');
        
        % Subplot for Investment
        subplot(1, 3, 2);
        hold on;
        title('Investment Policy');
        if fig == 1
            xlabel('Risky Asset');
        else
            xlabel('Safe Asset');
        end
        ylabel('Investment');
        
        % Subplot for Consumption
        subplot(1, 3, 3);
        hold on;
        title('Consumption Policy');
        if fig == 1
            xlabel('Risky Asset');
        else
            xlabel('Safe Asset');
        end
        ylabel('Consumption');
        
        legend_labels = cell(par.Nz, 1);
        
        for z = 1:par.Nz
            % Initialize arrays to store mean policy functions
            mean_sav = zeros(par.Na, par.Na_risk);
            mean_inv = zeros(par.Na, par.Na_risk);
            mean_cons = zeros(par.Na, par.Na_risk);
            
            % Calculate mean policy functions over all r
            for r = 1:par.Nr_risky
                mean_sav = mean_sav + reshape(sav_ret{z, r, t}, par.Na, par.Na_risk);
                mean_inv = mean_inv + reshape(inv_ret{z, r, t}, par.Na, par.Na_risk);
                mean_cons = mean_cons + reshape(cons_ret{z, r, t}, par.Na, par.Na_risk);
            end
            mean_sav = mean_sav / par.Nr_risky;
            mean_inv = mean_inv / par.Nr_risky;
            mean_cons = mean_cons / par.Nr_risky;
            
            if fig == 1 % Fixed median safe asset
                x_values = par.risk_agrid;
                sav_plot = mean_sav(median_safe_idx, :);
                inv_plot = mean_inv(median_safe_idx, :);
                cons_plot = mean_cons(median_safe_idx, :);
            else % Fixed median risky asset
                x_values = par.agrid;
                sav_plot = mean_sav(:, median_risk_idx);
                inv_plot = mean_inv(:, median_risk_idx);
                cons_plot = mean_cons(:, median_risk_idx);
            end
            
            % Plot Savings
            subplot(1, 3, 1);
            plot(x_values, sav_plot, 'Color', colors(z, :), 'LineWidth', 2);
            
            % Plot Investment
            subplot(1, 3, 2);
            plot(x_values, inv_plot, 'Color', colors(z, :), 'LineWidth', 2);
            
            % Plot Consumption
            subplot(1, 3, 3);
            plot(x_values, cons_plot, 'Color', colors(z, :), 'LineWidth', 2);
            
            % Create legend label
            legend_labels{z} = sprintf('z=%d', z);
        end
        
        % Add legend to the last subplot
        subplot(1, 3, 3);
        legend(legend_labels, 'Location', 'eastoutside');
        
        % Adjust subplot layouts to make room for the legend
        set(gcf, 'Position', [100, 100, 1500, 400]);
        
        % Add an overall title
        sgtitle(sprintf('Mean Policy Functions over r at t = %d, %s', par.Twork + t, figure_titles{fig}));
        
        % Save figure
        output_folder = 'Figures';
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        saveas(gcf, fullfile(output_folder, sprintf('mean_policy_functions_%s', strrep(figure_titles{fig}, ' ', '_'))));
        close(gcf);
    end
end
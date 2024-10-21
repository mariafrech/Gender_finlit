
% ------------------------------------------------------------------------
% Helper function to find optimal future asset allocation for working periods
function [v_max, idx] = max_future_value_work(budget, A_prime, z, r, t, par, V_w)
    % Calculate consumption for all possible future asset allocations
    cons = budget - A_prime(:,1) - A_prime(:,2);
    
    % Calculate utility
    u = par.U(cons);
    u(cons <= 0) = -inf;
    
    % Calculate expected future value
    EV = zeros(size(A_prime,1), 1);
    for zp = 1:par.Nz
        for rp = 1:par.Nr_risky
            EV = EV + par.Pi(z,zp) * par.Pi_risk(r,rp) * V_w{zp,rp,(t+1)};
        end
    end
    
    % Calculate total value
    v = u + par.beta * par.survprob(t) * EV;
    
    % Find maximum value and corresponding index
    [v_max, idx] = max(v);
end
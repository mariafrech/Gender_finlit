function [v_max, idx] = max_future_value_ret(budget, A_prime, z, r, t, par, V_ret)
    % Calculate consumption for all possible future asset allocations
    cons = budget - A_prime(:,1) - A_prime(:,2);                           %A_prime(:,1) safe assets (a1,...,a1,a2,...a2,..) and A_prime(:,2) risky assets (ar1,ar2,.....,ar1,...) 
    
    % Calculate utility
    u = par.U(cons);
    u(cons <= 0) = -inf;
    
    % Calculate expected future value
    EV = zeros(size(A_prime,1), 1);
    for rp = 1:par.Nr_risky
        EV = EV + par.Pi_risk(r,rp) * V_ret{z,rp,t+1};
    end

    

    

    % Calculate total value
    v = u + par.beta * par.survprob(par.Twork + t) * EV;
    
    % Find maximum value and corresponding index
    [v_max, idx] = max(v);
end



function CAIE()
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------

    clc;
    clear; close all;

    % Set default plot properties (SCI style: Times New Roman, thick lines)
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultLineLineWidth', 1.5);

    %% 1. System Parameter Initialization
    sys_params.WarrantyTime = 2000;    
    sys_params.TargetMTBF = 800;

    %  f_i 
    fi_vals = [0.25, 0.20, 0.20, 0.15, 0.25, 0.20];

    %  Boundaries
    subsystems(1) = struct('Name','Cutting(H)', 'Type',1, 'Beta',1.8, 'Tech_LB',2077.47, 'Tech_UB',3022.53, 'Abs_LB',0, 'Abs_UB',3089.32, 'Cm',1.20, 'Cd',1.55, 'Cr_Fixed',1.65,    'fi',fi_vals(1), 'ki_coeff',0.9);
    subsystems(2) = struct('Name','Seeding(H)', 'Type',1, 'Beta',1.6, 'Tech_LB',2375.57, 'Tech_UB',3624.43, 'Abs_LB',0, 'Abs_UB',3674.41, 'Cm',1.50, 'Cd',1.30, 'Cr_Fixed',1.60,    'fi',fi_vals(2), 'ki_coeff',0.8);
    subsystems(3) = struct('Name','Fertilizer(H)','Type',1,'Beta',2.8, 'Tech_LB',1559.37, 'Tech_UB',2340.63, 'Abs_LB',0, 'Abs_UB',2389.22, 'Cm',0.90, 'Cd',1.25, 'Cr_Fixed',1.55,    'fi',fi_vals(3), 'ki_coeff',0.85);
    subsystems(4) = struct('Name','Pressing(S)', 'Type',2, 'Beta',2.1, 'Tech_LB',1888.21, 'Tech_UB',2911.79, 'Abs_LB',0, 'Abs_UB',2988.21, 'Cm',2.30, 'Cd',1.55, 'Cr_Fixed',1.85, 'fi',fi_vals(4), 'ki_coeff',0.65);
    subsystems(5) = struct('Name','Drive(S)',    'Type',2, 'Beta',2.4, 'Tech_LB',2024.21, 'Tech_UB',3275.79, 'Abs_LB',0, 'Abs_UB',3358.76, 'Cm',2.10, 'Cd',1.50, 'Cr_Fixed',1.75, 'fi',fi_vals(5), 'ki_coeff',0.75);
    subsystems(6) = struct('Name','Frame(S)',    'Type',2, 'Beta',1.4, 'Tech_LB',1671.76, 'Tech_UB',2528.24, 'Abs_LB',0, 'Abs_UB',2521.34, 'Cm',3.05, 'Cd',1.70, 'Cr_Fixed',2.10, 'fi',fi_vals(6), 'ki_coeff',0.9);

    n = length(subsystems);
    lb = [subsystems.Tech_LB]; 
    ub = [subsystems.Tech_UB]; 

    %% 2. Particle Swarm Optimization (PSO)
    global history;
    history.fval = [];
    history.iteration = [];

    outputFcn = @(optimValues, state) recordHistory(optimValues, state);
    options = optimoptions('particleswarm', ...
        'SwarmSize', 200, ...         
        'MaxIterations', 300, ...    
        'Display', 'iter', ...        
        'OutputFcn', outputFcn);
        
    fun = @(theta) ObjectiveFunction(theta, subsystems, sys_params);
    
    fprintf('Running optimization...\n');
    [theta_opt, fval_opt] = particleswarm(fun, n, lb, ub, options);

    %% 3. Figures
    [TotalCost_Opt, C_mfg_opt, C_war_opt, sys_mtbf_opt, ~] = CalculateDetails(theta_opt, subsystems, sys_params);

    % --- Figure 1: Convergence ---
    figure('Color','w', 'Position', [50, 100, 550, 400]);
    plot(history.iteration, history.fval, 'b-', 'LineWidth', 2);
    xlabel('Iteration Number', 'Interpreter', 'latex');
    ylabel('Total Cost (CNY 10k)', 'Interpreter', 'latex');
    title('Convergence of PSO Algorithm', 'Interpreter', 'latex');
    grid on;
    
    inset_ax = axes('Position',[.6 .6 .25 .25]);
    box on;
    start_idx = max(1, length(history.fval)-50);
    plot(history.iteration(start_idx:end), history.fval(start_idx:end), 'b-');
    title('Final 50 Iterations');
    grid on;

    fprintf('\n======================================================\n');
    fprintf('>> Optimization completed!\n');
    fprintf('Final System MTBF: %.2f h (Target: %.0f h)\n', sys_mtbf_opt, sys_params.TargetMTBF);
    fprintf('Total Life Cycle Cost (LCC): %.4f (CNY 10k)\n', TotalCost_Opt);
    fprintf('======================================================\n');
end

%% Auxiliary function: record history
function stop = recordHistory(optimValues, state)
    stop = false;
    global history;
    switch state
        case 'iter'
            history.fval = [history.fval; optimValues.bestfval];
            history.iteration = [history.iteration; optimValues.iteration];
    end
end

%% Auxiliary function: objective function (including penalty)
function cost_penalty = ObjectiveFunction(theta, subs, params)
    [TotalCost, ~, ~, sys_mtbf, ~] = CalculateDetails(theta, subs, params);
    if sys_mtbf < params.TargetMTBF
        penalty = 1e6 * (params.TargetMTBF - sys_mtbf)^2;
    else
        penalty = 0;
    end
    cost_penalty = TotalCost + penalty;
end

%% Core calculation function (applying the latest economic allocation principles)
function [TotalCost, C_mfg_vec, C_war_vec, sys_mtbf, Failures_vec] = CalculateDetails(theta, subs, params)
    n = length(theta);
    C_mfg_vec = zeros(1, n); 
    C_war_vec = zeros(1, n);
    Failures_vec = zeros(1, n);
    T = params.WarrantyTime;
    
    for i = 1:n
        % --- 1. Manufacturing Cost ---
        A_Sigma = subs(i).Cm + subs(i).Cd;
        denominator = subs(i).Abs_UB - theta(i);
        if denominator <= 1e-4, denominator = 1e-4; end
        
        numerator = theta(i) - subs(i).Tech_LB;
        if numerator < 0, numerator = 0; end
        
        exponent = subs(i).fi * (numerator / denominator);
        C_mfg_vec(i) = A_Sigma * (1 + subs(i).ki_coeff * (exp(exponent) - 1));
        
        % --- 2. Warranty Cost Prediction and Calculation ---
        beta = subs(i).Beta;
        eta = theta(i) / gamma(1 + 1/beta);
        
        if subs(i).Type == 2 
            % Repairable subsystem
            ExpFailures = (T / eta)^beta;
            % Unit cost per repair.
            UnitCost = subs(i).Cr_Fixed;    
            
        else 
            % Replaceable subsystem (Renewal Process)
            % Subsystem is completely replaced upon failure.
            mu = theta(i);
            sigma = mu * sqrt(gamma(1+2/beta)/gamma(1+1/beta)^2 - 1);
            rho_i = 0.5 - (sigma/mu)^2 / 2;
            
            F_func = @(t) 1 - exp(-(t/eta)^beta);
            f_func = @(t) (beta/eta) * (t/eta)^(beta-1) * exp(-(t/eta)^beta);
            
            g_func = @(t) (F_func(t) / sqrt(mu * f_func(t) + 1e-8)) - t/mu + rho_i;
            
            options_fzero = optimset('Display', 'off');
            try
                ts = fzero(g_func, mu, options_fzero);
                if ts <= 0 || isnan(ts)
                    ts = mu; 
                end
            catch
                ts = mu; 
            end
            
            if T < ts
                F_ts = F_func(ts);
                denom_alpha = F_ts * (ts/mu - rho_i);
                if denom_alpha < 1e-6
                    denom_alpha = 1e-6; 
                end
                alpha_i = (ts/mu - rho_i - F_ts) / denom_alpha;
                
                F_T = F_func(T);
                ExpFailures = F_T / (1 - alpha_i * F_T);
            else
                ExpFailures = T / mu - rho_i;
            end
            
            % Unit cost for a replacement.
            UnitCost = subs(i).Cm + subs(i).Cd + subs(i).Cr_Fixed; 
        end
        
        Failures_vec(i) = ExpFailures;
        C_war_vec(i) = UnitCost * ExpFailures;
    end
    
    TotalCost = sum(C_mfg_vec) + sum(C_war_vec);
    
    % --- 3. System MTBF ---
    R_sys_func = @(t) 1;
    for k = 1:n
        eta_k = theta(k) / gamma(1 + 1/subs(k).Beta);
        R_sys_func = @(t) R_sys_func(t) .* exp( -(t./eta_k).^subs(k).Beta );
    end
    t_limit = params.TargetMTBF * 5;
    sys_mtbf = integral(R_sys_func, 0, t_limit);
end

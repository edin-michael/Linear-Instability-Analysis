%% Muller's method - Root finding algorithm with Mode Detection (Robust Iterative G Solver)
clear variables; close all; clc;

% --- 0. TARGET INPUT ---
fprintf('======================================\n');
fprintf('   Inverse Solver for Parameter G     \n');
fprintf('======================================\n');
target_BL = input('Enter the target Breakup Length (BL) in mm: ');
if isempty(target_BL) || target_BL <= 0
    error('Please enter a valid positive number for BL.');
end

% --- 1. Iteration Setup ---
% Initial guess for G
G = 1.07E-07; 

% Optimization parameters
max_G_iter = 50;       % Max attempts to find G
bl_tolerance = 0.001;    % Acceptable error in mm
G_history = [];        % To store history for Secant method
BL_history = [];       % To store history of results

% Flags for the outer loop
G_found = false;

fprintf('\nStarting iteration to find G for Target BL = %.2f mm...\n', target_BL);

% --- OUTER LOOP: Adjust G ---
for g_iter = 1:max_G_iter
    
    % --- 1. Physical Parameters (Recalculated with current G) ---
    m = 0; h1 = 1;
    % G is defined by the loop
    K_const = 0.36515; n = 0.453; sigma = 0.06875; rho_l = 1019.56;
    Qi = 1.188/rho_l; Qo = Qi; h = 0.576648413;
    We_l = 2665.909634; We_sl = 590.6091885;
    b = 0.001;
    We_i = 0; We_si = 0;
    We_o = 0; We_so = 0;

    U_l = sqrt(We_l*sigma/(rho_l*b));
    Re_n = rho_l*U_l^(2-n)*b^n/K_const;

    % --- 2. Solver Setup ---
    kvec = 0.01:0.01:25;
    Nk = numel(kvec);

    % Preallocate
    K_vec = NaN(Nk,1);
    S_sin_real = NaN(Nk,1);
    S_sin_imag = NaN(Nk,1);
    S_var_real = NaN(Nk,1);
    S_var_imag = NaN(Nk,1);

    % --- CONTINUATION VARIABLES INITIALIZATION ---
    S_seed = [0.0012 - 0.0040i, 0.0191 - 0.0650i, 0.0386 - 0.1320i];
    S_next_guess = S_seed; 
    low_growth_flag = false; 
    unstable_tol = 1e-10;
    sqrt_complex = @(z) sqrt(z);

    % --- 3. Loop through wavenumbers (Core Physics) ---
    for ik = 1:Nk
        k = kvec(ik);

        % Auxiliary Functions
        L = @(S) k.*sqrt_complex( Re_n.*(S + k.*1i) ./ (Re_n.*(S + k.*1i) + 2.*n.*k.^2*(2*G)^(n-1)) );
        M = @(S) k.*sqrt_complex( 1 + (4.*We_si./(We_l.*Qi)) ./ (S + k.*1i.*sqrt_complex(We_i./(We_l.*Qi))).^2 );

        % Boundary Condition Matrix Terms
        D1 = ((1-h^2.*k^2)./h^2 - We_sl./h^3 + h.*We_si)./We_l;
        D2 = ((1-k^2)./h^2 + We_sl + We_so)./We_l;

        D3 = @(S) (S + k.*1i).^2 .* ( K1L(m,L(S),h1).*I0Lh(m,L(S),h) + K0Lh(m,L(S),h).*I1L(m,L(S),h1) ) ...
                  ./ ( L(S) .* ( I1L(m,L(S),h1).*K1Lh(m,L(S),h) - K1L(m,L(S),h1).*I1Lh(m,L(S),h) ) );

        D4 = @(S) (S + k.*1i).^2 .* ( K1L(m,L(S),h1).*I0L(m,L(S),h1) + K0L(m,L(S),h1).*I1L(m,L(S),h1) ) ...
                  ./ ( L(S) .* ( I1L(m,L(S),h1).*K1Lh(m,L(S),h) - K1L(m,L(S),h1).*I1Lh(m,L(S),h) ) );

        D5 = @(S) (S + k.*1i).^2 .* ( K1Lh(m,L(S),h).*I0Lh(m,L(S),h) + K0Lh(m,L(S),h).*I1Lh(m,L(S),h) ) ...
                  ./ ( L(S) .* ( I1L(m,L(S),h1).*K1Lh(m,L(S),h) - K1L(m,L(S),h1).*I1Lh(m,L(S),h) ) );

        D6 = @(S) (S + k.*1i).^2 .* ( K1Lh(m,L(S),h).*I0L(m,L(S),h1) + K0L(m,L(S),h1).*I1Lh(m,L(S),h) ) ...
                  ./ ( L(S) .* ( I1L(m,L(S),h1).*K1Lh(m,L(S),h) - K1L(m,L(S),h1).*I1Lh(m,L(S),h) ) );

        D7 = @(S) (1./M(S)) .* ( (S.*Qi.^0.5 + 1i.*k.*sqrt_complex(We_i./We_l)).^2 + 4.*We_si./We_l ) .* ( I0Mh(m,M(S),h)./I1Mh(m,M(S),h) );
        D8 = @(S) (1./k) .* (S.*Qo.^0.5 + 1i.*k.*sqrt_complex(We_o./We_l)).^2 .* ( K0k(m,k,h1)./K1k(m,k,h1) );

        % Characteristic function
        f = @(S) ( D3(S) + D7(S) - D1 ) .* ( D6(S) + D8(S) - D2 ) - D4(S).*D5(S);

        % Divided differences
        f2 = @(S1,S2) ( f(S1) - f(S2) ) ./ (S1 - S2);
        f3 = @(S1,S2,S3) ( f2(S1,S2) - f2(S2,S3) ) ./ (S1 - S3);

        % Muller's Method Initialization
        Nmax = 1000; tol = 1e-8;
        S_old = S_next_guess;
        S_root = NaN; converged = false;
        
        for iter = 1:Nmax
            try
                Om = f2(S_old(3), S_old(2)) + f2(S_old(3), S_old(1)) - f2(S_old(2), S_old(1));
                discr = Om^2 - 4 * f(S_old(3)) * f3(S_old(3), S_old(2), S_old(1));
                Discrim = sqrt_complex(discr);
                A1 = Om + Discrim; A2 = Om - Discrim;
                if abs(A1) >= abs(A2), A = A1; else, A = A2; end

                S_new = S_old(3) - 2*f(S_old(3)) / A;
                if ~isfinite(S_new) || isnan(real(S_new)) || isnan(imag(S_new)), break; end
                if abs(S_new - S_old(3)) < tol, S_root = S_new; converged = true; break; end
                S_old = [S_old(2), S_old(3), S_new];
            catch ME, break; end
        end

        if ~converged
            if exist('S_new','var') && isfinite(S_new)
                S_root = S_new;
            else
                if ~low_growth_flag, S_next_guess = S_seed; end
                continue; 
            end
        end

        Sreal = real(S_root); Simag = imag(S_root);

        if Sreal > unstable_tol
            S_new_guess = S_root; 
            if Sreal <= 0.001, low_growth_flag = true; else, low_growth_flag = false; end
            S_next_guess = [S_next_guess(2), S_next_guess(3), S_new_guess];

            try
                Val_Numerator = D3(S_root) + D7(S_root) - D1;
                Val_Denominator = D5(S_root);
                eta_ratio = Val_Numerator ./ Val_Denominator;
            catch, eta_ratio = NaN; end

            K_vec(ik) = k;
            if ~isnan(real(eta_ratio)) && real(eta_ratio) > 0
                S_sin_real(ik) = Sreal; S_sin_imag(ik) = Simag;
            else
                S_var_real(ik) = Sreal; S_var_imag(ik) = Simag;
            end
        else
            if low_growth_flag, S_next_guess = S_next_guess; else, S_next_guess = [S_next_guess(2), S_next_guess(3), S_root]; end
        end
        clear S_new S_root S_old Om Discrim A1 A2 A Val_Numerator Val_Denominator eta_ratio
    end

    % --- 4. CALCULATE RESULTING BL ---
    All_Sreal = [S_sin_real; S_var_real];
    All_Sreal = All_Sreal(~isnan(All_Sreal));
    
    if isempty(All_Sreal)
        current_BL = NaN;
        Sr_max = 0;
    else
        [Sr_max, ~] = max([S_sin_real; S_var_real]);
        sr_max_dim = Sr_max * U_l / b;
        if sr_max_dim ~= 0
            current_BL = U_l / sr_max_dim * 3 * 1000; % mm
        else
            current_BL = NaN;
        end
    end
    
    % --- 5. CHECK CONVERGENCE & UPDATE G ---
    fprintf('Iter %d: G = %.4e -> BL = %.4f mm (Error: %.4f)\n', ...
        g_iter, G, current_BL, current_BL - target_BL);
    
    if isnan(current_BL)
        fprintf('  > Solver failed to find instability. Increasing G slightly.\n');
        G = G * 1.5; % Boost G to find instability
        continue;
    end
    
    % Check tolerance
    if abs(current_BL - target_BL) < bl_tolerance
        fprintf('>> CONVERGENCE ACHIEVED!\n');
        G_found = true;
        break; 
    end
    
    % Store history
    G_history(end+1) = G;
    BL_history(end+1) = current_BL;
    
    % Update G using Secant Method logic
    if length(G_history) == 1
        % Only one point, we need a second point to calculate slope.
        % We know increasing G decreases BL.
        if current_BL > target_BL
            % BL is too big -> We need smaller BL -> Increase G
            G = G * 1.2; 
        else
            % BL is too small -> We need larger BL -> Decrease G
            G = G * 0.8;
        end
    else
        % We have at least 2 points. Use Secant Method.
        % G_new = G_curr - f(G_curr) * (G_curr - G_prev) / (f(G_curr) - f(G_prev))
        % where f(G) = BL_calc - BL_target
        
        G_curr = G_history(end);
        G_prev = G_history(end-1);
        BL_curr = BL_history(end);
        BL_prev = BL_history(end-1);
        
        denom = (BL_curr - BL_prev);
        
        if abs(denom) < 1e-6
            % Avoid division by zero if BL didn't change
            G = G * 1.05; 
        else
            G_new = G_curr - (BL_curr - target_BL) * (G_curr - G_prev) / denom;
            
            % Damping/Bounds Check to prevent negative G or wild swings
            if G_new <= 0
                G_new = G_curr * 0.5; % Just halve it if it tries to go negative
            end
            
            G = G_new;
        end
    end
end

if ~G_found
    warning('Maximum iterations reached. Showing results for last calculated G.');
end

% --- 6. Post processing & plotting (Uses data from final iteration) ---
has_sinuous = any(~isnan(S_sin_real));
has_varicose = any(~isnan(S_var_real));

fprintf('\n======================================\n');
fprintf('FINAL RESULTS FOR G = %.4e\n', G);
if has_sinuous, disp('>> Para-Sinuous Mode: PRESENT'); else disp('>> Para-Sinuous Mode: NOT PRESENT'); end
if has_varicose, disp('>> Para-Varicose Mode: PRESENT'); else disp('>> Para-Varicose Mode: NOT PRESENT'); end
fprintf('--------------------------------------\n');

if ~isempty(All_Sreal)
    [Sr_max, max_idx] = max([S_sin_real; S_var_real]);
    
    if max_idx <= Nk 
        k_max_nd = K_vec(max_idx); Simag_max_nd = S_sin_imag(max_idx);
    else 
        idx_var = max_idx - Nk;
        k_max_nd = K_vec(idx_var); Simag_max_nd = S_var_imag(idx_var);
    end
    
    k_dim_max = k_max_nd / b; 
    lambda_max = 2 * pi / k_dim_max * 1000000;
    omega_i_dim_max = Simag_max_nd * U_l / b; 
    f_dist_max = omega_i_dim_max / (2 * pi); 
    Vp_max = omega_i_dim_max / k_dim_max; 
    
    % Plot
    figure(1); clf; hold on; box on; grid on;
    if has_sinuous
        plot(K_vec, S_sin_real, 'b-o', 'LineWidth',1,'MarkerSize',4,'DisplayName','Para-Sinuous');
    end
    if has_varicose
        plot(K_vec, S_var_real, 'r--s', 'LineWidth',1,'MarkerSize',4,'DisplayName','Para-Varicose');
    end
    xlabel('Wavenumber (k)');
    ylabel('Growth Rate (S_{real})');
    title(sprintf('Dispersion Curve (G=%.2e, BL=%.2fmm)', G, current_BL));
    legend('show','Location','Best');

    fprintf('Maximum Growth Rate (non-dim): %.4f\n', Sr_max);
    fprintf('Est. Breakup Length Parameter (BL): %.2f (mm)\n', current_BL);
    fprintf('Dominant Wavelength (lambda): %.2f (mum)\n', lambda_max);
    fprintf('Disturbance Frequency (f_dist): %.4f (Hz)\n', f_dist_max);
    fprintf('======================================\n');
else
    disp('No unstable modes found for this G.');
end
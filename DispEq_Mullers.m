%% Muller's method - Root finding algorithm with Mode Detection (robust version)
clear variables; close all; clc;

% --- 1. Physical Parameters ---
m = 0; h1 = 1;
G = 1.88E-6;
K_const = 0.000907919; n = 1; sigma = 0.072; rho_l = 997.38;
Qi = 1.188/rho_l; Qo = Qi; h = 0.581247993;
d = 0.002; % orifice dia, m
b = 0.5 * d; % orifice radius, m
% m_dot = 0.0056; % mass flow rate, kg/s
% V_dot = m_dot/rho_l; % volume flow rate
U_l = 13.65013612; %V_dot/(pi/4*d^2);

We_l = rho_l * U_l^2 * b / sigma; 
We_sl = 637.3443705;
We_i = 0; We_si = 0;
We_o = 0; We_so = 0;

% U_l = sqrt(We_l*sigma/(rho_l*b));
Re_n = rho_l*U_l^(2-n)*b^n/K_const;

% --- 2. Solver Setup ---
kvec = 0.01:0.01:50;
Nk = numel(kvec);

% Preallocate (store only unstable entries; max possible = Nk)
K_vec = NaN(Nk,1);
S_sin_real = NaN(Nk,1);
S_sin_imag = NaN(Nk,1);
S_var_real = NaN(Nk,1);
S_var_imag = NaN(Nk,1);

% --- CONTINUATION VARIABLES INITIALIZATION ---
% 1. Initial guess seed values (only used for the very first iteration)
S_seed = [0.0012 - 0.0040i, 0.0191 - 0.0650i, 0.0386 - 0.1320i];
S_next_guess = S_seed; % The current set of guesses for S_old
% 2. Threshold for low growth regime where aggressive search is stopped
low_growth_flag = false; 
% ---------------------------------------------

% small tolerance for calling an S "unstable"
unstable_tol = 1e-10;

% Helper: explicit sqrt wrapper (keeps code readable; uses MATLAB sqrt which handles complex)
sqrt_complex = @(z) sqrt(z);

% --- 3. Loop through wavenumbers ---
for ik = 1:Nk
    k = kvec(ik);

    % --- Auxiliary Functions (take care: L and M will be complex functions of S) ---
    L = @(S) k.*sqrt_complex( Re_n.*(S + k.*1i) ./ (Re_n.*(S + k.*1i) + 2.*n.*k.^2*(2*G)^(n-1)) );
    M = @(S) k.*sqrt_complex( 1 + (4.*We_si./(We_l.*Qi)) ./ (S + k.*1i.*sqrt_complex(We_i./(We_l.*Qi))).^2 );

    % Boundary Condition Matrix Terms (note: these call external Bessel wrappers)
    D1 = ((1-h^2.*k^2)./h^2 - We_sl./h^3 + h.*We_si)./We_l;
    D2 = ((1-k^2)./h^2 + We_sl + We_so)./We_l;

    % The following D3..D6 and D7..D8 are function handles of S.
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

    % Divided differences (first and second)
    f2 = @(S1,S2) ( f(S1) - f(S2) ) ./ (S1 - S2); % clearer form
    f3 = @(S1,S2,S3) ( f2(S1,S2) - f2(S2,S3) ) ./ (S1 - S3);

    % --- Muller's Method Initialization ---
    Nmax = 1000;
    tol = 1e-8;

    % Use the guesses from the previous iteration (continuation)
    S_old = S_next_guess;
    
    S_root = NaN;
    converged = false;
    
    % iterate
    for iter = 1:Nmax
        try
            Om = f2(S_old(3), S_old(2)) + f2(S_old(3), S_old(1)) - f2(S_old(2), S_old(1));
            discr = Om^2 - 4 * f(S_old(3)) * f3(S_old(3), S_old(2), S_old(1));
            Discrim = sqrt_complex(discr);

            A1 = Om + Discrim;
            A2 = Om - Discrim;
            if abs(A1) >= abs(A2)
                A = A1;
            else
                A = A2;
            end

            S_new = S_old(3) - 2*f(S_old(3)) / A;
            if ~isfinite(S_new) || isnan(real(S_new)) || isnan(imag(S_new))
                break; % fail this k
            end

            if abs(S_new - S_old(3)) < tol
                S_root = S_new;
                converged = true;
                break;
            end

            % shift for next iteration
            S_old = [S_old(2), S_old(3), S_new];

        catch ME
            % numerical failure (e.g., evaluation of Bessel wrappers) -> break and mark fail
            % warning('k=%.4f: Mullers step failed at iter %d: %s', k, iter, ME.message);
            break;
        end
    end

    if ~converged
        % Fallback for non-convergence
        if exist('S_new','var') && isfinite(S_new)
            S_root = S_new;
            % warning('k=%.4f: Mullers did not converge to tol but produced a finite root (iter=%d).', k, iter);
        else
            % If the solver fails completely AND we are NOT in the low-growth phase,
            % try re-seeding the guess from the original S_seed for better stability.
            if ~low_growth_flag
                 % OPTIONAL: Add aggressive perturbation search here if simple continuation fails.
                 S_next_guess = S_seed; % Reset guess to the original seed
            end
            
            continue; % skip to next k; S_root remains NaN.
        end
    end

    Sreal = real(S_root);
    Simag = imag(S_root);

    % --- Mode classification --- 
    if Sreal > unstable_tol
        
        % --- UPDATE CONTINUATION GUESSES ---
        % Only update the continuation guess if a successful unstable root is found
        S_new_guess = S_root; 
        
        % Check the low growth threshold
        if Sreal <= 0.001
            low_growth_flag = true;
        else
            low_growth_flag = false;
        end
        
        % Prepare guess for next iteration (ik+1)
        S_next_guess = [S_next_guess(2), S_next_guess(3), S_new_guess];


        % --- CLASSIFICATION (CORRECTED RATIO) ---
        % compute amplitude ratio eta_out/eta_in = (D3+D7-D1)/D5
        try
            Val_Numerator = D3(S_root) + D7(S_root) - D1;
            Val_Denominator = D5(S_root);
            eta_ratio = Val_Numerator ./ Val_Denominator;
        catch
            eta_ratio = NaN;
        end

        % Store at the current k index (ik)
        K_vec(ik) = k;

        % classify by real(eta_ratio) sign
        if ~isnan(real(eta_ratio)) && real(eta_ratio) > 0
            S_sin_real(ik) = Sreal;
            S_sin_imag(ik) = Simag;
        else
            S_var_real(ik) = Sreal;
            S_var_imag(ik) = Simag;
        end

        fprintf('k=%.4f : Unstable root found (S=%.4f %+.2fi)\n', k, Sreal, Simag);
    else
        % If the root is stable, we stop propagating the guess if we hit the low-growth threshold
        if low_growth_flag
            S_next_guess = S_next_guess; % Keep the last guess for stability
        else
            % If it's stable but NOT past the low-growth point,
            % use the stable root as the guess for the next iteration anyway.
            S_next_guess = [S_next_guess(2), S_next_guess(3), S_root];
        end
    end

    % clear iteration-specific vars to avoid accidental reuse
    clear S_new S_root S_old Om Discrim Discrim A1 A2 A Val_Numerator Val_Denominator eta_ratio
end

% --- Post processing & plotting ---
has_sinuous = any(~isnan(S_sin_real));
has_varicose = any(~isnan(S_var_real));

fprintf('\n--------------------------------------\n');
if has_sinuous, disp('>> Para-Sinuous Mode: PRESENT'); else disp('>> Para-Sinuous Mode: NOT PRESENT'); end
if has_varicose, disp('>> Para-Varicose Mode: PRESENT'); else disp('>> Para-Varicose Mode: NOT PRESENT'); end
fprintf('--------------------------------------\n');

All_Sreal = [S_sin_real; S_var_real];
All_Sreal = All_Sreal(~isnan(All_Sreal));

if isempty(All_Sreal)
    disp('No unstable modes found.');
else
    % --- 6. FIND MAXIMUM GROWTH RATE MODE ---
    
    % Combine all unstable growth rates and find the index of the maximum
    [Sr_max, max_idx] = max([S_sin_real; S_var_real]);
    
    % Determine which array the maximum came from and get the corresponding k and Simag
    if max_idx <= Nk % Max came from S_sin_real
        k_max_nd = K_vec(max_idx);
        Simag_max_nd = S_sin_imag(max_idx);
    else % Max came from S_var_real
        % Adjust index to the original Nk-sized array (since S_var_real starts after S_sin_real)
        idx_var = max_idx - Nk;
        k_max_nd = K_vec(idx_var);
        Simag_max_nd = S_var_imag(idx_var);
    end
    
    % Handle potential NaN/finite issue if max_idx points to NaN (shouldn't happen, but safety first)
    if isnan(Sr_max) || isnan(k_max_nd)
         disp('Error finding corresponding max k/Simag. Check storage vectors.');
         return;
    end
    
    % --- 7. DIMENSIONAL CALCULATIONS FOR MAX MODE ---
    
    % A. Dimensional Wavenumber
    k_dim_max = k_max_nd / b; 
    
    % B. Wavelength (lambda)
    lambda_max = 2 * pi / k_dim_max * 1000000;
    
    % C. Dimensional Imaginary Frequency (omega_i, dim)
    omega_i_dim_max = Simag_max_nd * U_l / b; 
    
    % D. Disturbance Frequency (f_dist)
    f_dist_max = omega_i_dim_max / (2 * pi); 
    
    % E. Propagation Velocity (Vp)
    Vp_max = omega_i_dim_max / k_dim_max; % Equivalent to (Simag_max_nd / k_max_nd) * U_l
    
    % --- 8. OUTPUT AND PLOTTING ---
    
    % Plot (as in your original code)
    figure(1); clf; hold on; box on; grid on;
    if has_sinuous
        plot(K_vec, S_sin_real, 'b-o', 'LineWidth',1,'MarkerSize',4,'DisplayName','Para-Sinuous');
    end
    if has_varicose
        plot(K_vec, S_var_real, 'r--s', 'LineWidth',1,'MarkerSize',4,'DisplayName','Para-Varicose');
    end
    xlabel('Wavenumber (k)');
    ylabel('Growth Rate (S_{real})');
    title(sprintf('Dispersion Curve (Max Growth Rate = %.4e)', Sr_max));
    legend('show','Location','Best');

    sr_max_dim = Sr_max * U_l / b;
    
    % note: ensure Sr_max not zero
    if sr_max_dim ~= 0
        BL = U_l / sr_max_dim * 12 * 1000; % keep your original formula (units check required)
    else
        BL = NaN;
    end

    fprintf('\n======================================\n');
    fprintf('DIMENSIONAL RESULTS (For Max Growth Rate Mode)\n');
    fprintf('--------------------------------------\n');
    fprintf('Maximum Growth Rate (non-dim): %.4f\n', Sr_max);
    fprintf('Maximum Growth Rate (dim) = Sr_max*U_l/b: %.4f (s^-1)\n', sr_max_dim);
    fprintf('Est. Breakup Length Parameter (BL): %.2f (mm)\n', BL);
    fprintf('--------------------------------------\n');
    fprintf('Dominant Wavenumber (k_dim): %.4f (m^-1)\n', k_dim_max);
    fprintf('Dominant Wavelength (lambda): %.2f (mum)\n', lambda_max);
    fprintf('Disturbance Frequency (f_dist): %.4f (Hz)\n', f_dist_max);
    fprintf('Propagation Velocity (Vp): %.4f (m/s)\n', Vp_max);
    fprintf('======================================\n');
end


%% ===== SAVE RESULTS TO EXCEL (Choose location manually) =====

% Ask user to choose location and file name
[FileName, PathName] = uiputfile('*.xlsx', 'Save Results As');
if isequal(FileName,0)
    disp('User cancelled file save.');
else
    SaveFile = fullfile(PathName, FileName);

    % Prepare table for vector data
    T = table(K_vec, S_sin_real, S_sin_imag, ...
              'VariableNames', {'K_vec', 'S_sin_real', 'S_sin_imag'});

    % Write vector data to Excel
    writetable(T, SaveFile, 'Sheet', 'Data', 'Range', 'A1');

    % Write scalar summary values
    summaryNames = {'BL'; 'lambda_max'; 'f_dist_max'; 'Vp_max'};
    summaryValues = [BL; lambda_max; f_dist_max; Vp_max];

    SummaryTable = table(summaryNames, summaryValues, ...
        'VariableNames', {'Parameter', 'Value'});

    writetable(SummaryTable, SaveFile, 'Sheet', 'Summary', 'Range', 'A1');

    fprintf('Results successfully saved to:\n%s\n', SaveFile);
end

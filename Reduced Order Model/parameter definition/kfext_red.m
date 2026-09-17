%% Kalman aumentado con estado r para modelo reducido
% r = (V_bus / V0)²  →  r=1 nominal, r<1 cuando cae el bus

% Una sola columna de acoplamiento (discreta, modelo reducido)
col_r_red = Bd_u_red * Ufeedforward_red;  % (8×1)

% Sistema aumentado (9×9)
Ad_aug_kf_red = [Ad_red,          col_r_red; ...
                 zeros(1,N_red),  1        ];   % r[k+1] = r[k]

C_aug_kf_red = [Cd_red, zeros(2,1)];            % (2×9)
B_aug_kf_red = [Bd_u_red; zeros(1,2)];          % (9×2)

% Ruido de proceso para r
dVdt_real = 4/60;
V_inicial  = 34;
dr_dt_max  = 2 * V_inicial / 34^2 * dVdt_real;
sigma_r    = 0.01;
fprintf('sigma_r: %.6f\n', sigma_r)

% Matrices de ruido aumentadas (Qn_red ya es discreto: Qnc_red*Ts)
Qn_ext_red = blkdiag(Qn_red, sigma_r^2);  % (9×9)
Rn_ext_red = Rn_red;                       % (2×2)

% Kalman discreto (dare)
[P_ext_red, ~, ~] = dare(Ad_aug_kf_red', C_aug_kf_red', Qn_ext_red, Rn_ext_red);
Kf_ext_red = (P_ext_red * C_aug_kf_red') / ...
             (C_aug_kf_red * P_ext_red * C_aug_kf_red' + Rn_ext_red);

% Matrices para bloque Simulink (State-Space discreto)
% Entradas: [u1, u2, y1, y2]  (4 entradas)
% Salidas:  x̂ (9×1) — estados 1:8 = estados térmicos, estado 9 = r̂
A_kf_red = Ad_aug_kf_red - Kf_ext_red * C_aug_kf_red;  % (9×9)
B_kf_red = [B_aug_kf_red, Kf_ext_red];                 % (9×4)
C_kf_red = eye(N_red+1);                                % (9×9)
D_kf_red = zeros(N_red+1, 4);                           % (9×4)

syskf_ext_red= ss(Ad_aug_kf_red, B_aug_kf_red, C_aug_kf_red, D_ss_red, Ts);

% Condición inicial: Δx=0 (desviaciones), r=1 (voltaje nominal)
x0_kf_red = [zeros(N_red,1); 1];   % (9×1)

% Verificaciones
fprintf('Rango observabilidad: %d / %d\n', ...
        rank(obsv(Ad_aug_kf_red, C_aug_kf_red)), N_red+1)
fprintf('Kf_ext_red fila r̂: [%.6f, %.6f]\n', ...
        Kf_ext_red(N_red+1,1), Kf_ext_red(N_red+1,2))
fprintf('|polo|_max estimador: %.8f\n', max(abs(eig(A_kf_red))))

% Verificar estabilidad
if max(abs(eig(A_kf_red))) < 1
    fprintf('✓ Estimador estable\n')
else
    fprintf('✗ Estimador inestable — revisar Qn_ext_red / Rn_ext_red\n')
end

% Mostrar ganancia de r̂ para verificar que es observable
if abs(Kf_ext_red(N_red+1,1)) < 1e-8 && abs(Kf_ext_red(N_red+1,2)) < 1e-8
    fprintf('⚠ Ganancia Kf para r̂ ≈ 0 — r puede no ser observable\n')
    fprintf('  Considera aumentar sigma_r\n')
else
    fprintf('✓ Ganancia Kf r̂: [%.8f, %.8f]\n', ...
            Kf_ext_red(N_red+1,1), Kf_ext_red(N_red+1,2))
end

dVdt = -0;
Vb0 = 30;
V0 = 30;
steptime = 20;
step_temp = 1;
ramp_wd = 0;

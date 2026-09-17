%% para aumentar con r filtro de kalman (NO es lo mismo que Kalman extendido)

% r escalar
% r = (V_bus_real / V_bus_nom)²

% Una sola columna de acoplamiento (discreta)
col_r = Bd_u * Ufeedforward;  % (38×1)  ← Bd_u en vez de B_u

% Sistema aumentado (39×39)
Ad_aug_kf = [Ad,          col_r;
             zeros(1,N),  1    ];   % r[k+1] = r[k]  ← 1 en vez de 0

C_aug_kf = [Cd, zeros(2,1)];       % (2×39)  ← Cd en vez de C_ss

B_aug_kf = [Bd_u;
            zeros(1,2)];            % (39×2)  ← Bd_u en vez de B_u

% Ruido de proceso

sigma_r = 0.01;
fprintf('sigma_r: %.6f\n', sigma_r)

Qn_ext = blkdiag(Qn, sigma_r^2);   % (39×39)  ← Qn ya es discreto (Qnc*Ts)
Rn_ext = Rn;                         % (2×2)   ← Rn ya es discreto (Rnc/Ts)

% Kalman discreto (dare en vez de lqe)
[P_ext, ~, ~] = dare(Ad_aug_kf', C_aug_kf', Qn_ext, Rn_ext);
Kf_ext = (P_ext * C_aug_kf') / (C_aug_kf * P_ext * C_aug_kf' + Rn_ext);

% Matrices bloque Simulink (igual que antes)
A_kf = Ad_aug_kf - Kf_ext * C_aug_kf;   % (39×39)
B_kf = [B_aug_kf, Kf_ext];              % (39×4)
C_kf = eye(N+1);                         % (39×39)
D_kf = zeros(N+1, 4);                   % (39×4)

syskf_ext= ss(Ad_aug_kf, B_aug_kf, C_aug_kf, D_ss, Ts);

dVdt = 0;
Vb0 = 30;
V0 = 30;
steptime = 20;
step_temp = 1;
wd = 0;

% Verificaciones
fprintf('Rango observabilidad: %d / %d\n', rank(obsv(Ad_aug_kf, C_aug_kf)), N+1)
fprintf('Kf_ext fila r̂: [%.6f, %.6f]\n', Kf_ext(N+1,1), Kf_ext(N+1,2))
fprintf('|polo|_max estimador: %.8f\n', max(abs(eig(A_kf))))


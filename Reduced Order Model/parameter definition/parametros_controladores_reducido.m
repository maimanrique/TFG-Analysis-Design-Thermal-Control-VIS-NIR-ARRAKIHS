%% 1. Cargar y procesar modelo
esatanmodel = readThermalModel_red("ARRAKIHS_ACASE2_reduced", "modelo_reducido_tfg/Modelo_reducido");
tmm_red = processModel_red(esatanmodel);

Q0_linear_red = tmm_red.Q0_linear;
Q0_red = tmm_red.Q0;
%% 2. Construir matrices de estado-espacio linealizado

N_red = tmm_red.N;  % Número de nodos D (difusión)
C_red = diag(tmm_red.C);  % Diag de capacidades
A_red = C_red \ tmm_red.K0;  % A = C**-1 * K0



%% 4. Sistema MIMO: u heater en Q(nodeIdx), y = T(nodeIdx)

% Índices de los nodos 
idxH1_red = 3;     % Heater 1 (radiador)
idxH2_red = 5;    % Heater 2 (cold finger)

% Matriz B (2 columnas: una por cada controlador)
B_red = zeros(N_red, 2); 

B_red(idxH1_red, 1) = 1 ./ tmm_red.C(idxH1_red);
B_red(idxH2_red, 2) = 1 ./ tmm_red.C(idxH2_red); 


%matriz identidad para comprobar x
CI_red = eye(8);
DI_red = zeros(8,3);
DII_red = zeros(8,2);

% Matriz C_ss (2 filas: una para cada sensor/medida)
C_ss_red = zeros(2, N_red); 

C_ss_red(1, idxH1_red) = 1;                     
C_ss_red(2, idxH2_red) = 1;

D_ss_red = zeros(2, 2);



T0_red = tmm_red.T0;  % Temperatura setpoint inicial = T0

%  B original  para los heaters (u)
B_u_red = zeros(N_red, 2);

B_u_red(idxH1_red, 1) = 1 ./ tmm_red.C(idxH1_red);  %  H1
B_u_red(idxH2_red, 2) = 1 ./ tmm_red.C(idxH2_red);  % H2 

% Creamos una B para las perturbaciones en el radiador(wd)
% Si wd entra en todos los nodos, Bd es la identidad dividida por la capacidad

idx_rad_red = 1:3;  
B_d_red = zeros(N_red,1);
B_d_red(idx_rad_red,1) = (1/24) ./ tmm_red.C(idx_rad_red);


% Ahora la entrada U del bloque será un vector [u; wd]
B_total_red = [B_u_red, B_d_red]; 


% La matriz D debe  coincidir con las dos entradas
D_total_red = [D_ss_red, zeros(size(C_ss_red, 1), size(B_d_red, 2))];

sys_red = ss(A_red, B_total_red, C_ss_red, D_total_red);


%% 5. PI 
Kp1 = 60 ;    Ki1 = 0.5;
Kp2 = 2 ;    Ki2 = 0.5;

%% parámetros
Pmin = 5; Vmin = 28;
V0 = 30; 
R0 = (Vmin^2)/Pmin; 
P0 = (V0^2)/R0;
Vmax = 34;
Pmax = (Vmax^2)/R0;
% 
% dVdt = -4/60; % de 34 a 28 +-4v/min 
dVdt = 0;
Vb0 = 30;


TH2_red = T0_red(idxH2_red);
TH_red = T0_red(idxH1_red);
Tsp_red = [TH_red; TH2_red];

RH1 = R0;
RH2 = R0;

Q0_lin_red = [Q0_linear_red(idxH1_red); Q0_linear_red(idxH2_red)];
Ufeedforward_red = [Q0_red(idxH1_red) ; Q0_red(idxH2_red)];

step_temp = 0; % 1 o 2 según caso ensayado
steptime = 20;

wd = -1;



%% ruido
sigma = 0.0025;


%% sistema aumentado

ny_red = size(C_ss_red,1);

% %continuo
% A_aug = [A_red zeros(N_red,ny_red);
%         -C_ss_red zeros(ny_red)];
% 
% B_aug = [B_u_red;
%          zeros(ny_red,2)];

%% discretización
Ts = 1;

sys_d_red = c2d(sys_red,   Ts, 'zoh');

[Ad_red, Bd_red, Cd_red, Dd_red] = ssdata(sys_d_red);  

sys_d_full_red = c2d(ss(A_red, [B_u_red, B_d_red], C_ss_red, zeros(ny_red,3)), Ts, 'zoh');

Bd_u_red  = sys_d_full_red.B(:,1:2);   % parte de heaters
Bd_d_red  = sys_d_full_red.B(:,3);     % parte de perturbación


Ad_aug_red = [Ad_red,          zeros(N_red, ny_red);
          -Cd_red*Ts,      eye(ny_red)     ];

Bd_aug_red = [Bd_u_red;
          zeros(ny_red,2) ];



%% LQR

% Q_lqr = diag([ones(1,N)*0.001, 1, 10]);
% R_lqr = diag([100, 50]);
Q_x_red = ones(1,N_red) * 0.001;
Q_x_red(idxH1_red) = 0.5;    
Q_x_red(idxH2_red) = 0.25;
Q_lqr_red = diag([Q_x_red, 50, 25]); 

% Q_lqr(idxH1, idxH1) = diag([1, 1]);          
% Q_lqr(idxH2, idxH2) = 1;                     

R_lqr_red = diag([3 , 12 ]);




%discretizado 
K_total_red = dlqr(Ad_aug_red, Bd_aug_red, Q_lqr_red, R_lqr_red);
K_red    = K_total_red(:, 1:N_red);       % ganancia de estado
Kint_red = K_total_red(:, N_red+1:N_red+ny_red);  % ganancia integral

% % continuo
% 
% K_total = lqr(A_aug, B_aug, Q_lqr_red, R_lqr_red);
% K = K_total(:,1:N_red);
% Kint = K_total(:,N_red+1:N_red+2);

%% Kalman

Qnc_red = 1e-3*eye(N_red);
Qnc_red(idx_rad_red,idx_rad_red)= 1e-2*eye(3);
Qnc_red(idxH1_red,idxH1_red) = 1e-5;
Qnc_red(idxH2_red,idxH2_red)    = 1e-5;

Rnc_red = 6.25e-6 * eye(2);

Qn_red = Qnc_red;
Rn_red = Rnc_red;


% % continuo
% Kf = lqe(A_red, eye(N_red),C_ss_red, Qn_red, Rn_red );
% syskf= ss(A_red, B_u_red, C_ss_red, D_ss_red);


% discretizado

Qn_red = Qnc_red*Ts;
Rn_red = Rnc_red/Ts;

[P_kf_d, ~, Kf_d_raw] = dare(Ad_red', Cd_red', Qn_red, Rn_red);

Kf = (P_kf_d * Cd_red') / (Cd_red * P_kf_d * Cd_red' + Rn_red);
syskf_red = ss(Ad_red, Bd_u_red, Cd_red, D_ss_red, Ts);

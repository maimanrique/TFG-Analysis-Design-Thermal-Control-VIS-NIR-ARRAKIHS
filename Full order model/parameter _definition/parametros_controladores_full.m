%comment o uncomment subsecciones en continuo o discreto en función del
%caso a estudiar


%% 1. Cargar y procesar modelo
esatanmodel = readThermalModel("ARRAKIHS_ACASE2", "Thermal_models/Simple_model");
tmm = processModel(esatanmodel);

Q0_linear = tmm.Q0_linear;
Q0 = tmm.Q0;

%% 2. Construir matrices de estado-espacio linealizado

N = tmm.N;  % Número de nodos D (difusión)
C = diag(tmm.C);  % Diag de capacidades
A = C \ tmm.K0;  % A = C**-1 * K0



%% 4. Sistema MIMO: u heater en Q(nodeIdx), y = T(nodeIdx)

% Índices de los nodos (ajusta según tu numeración en tmm)
idxH1 = [23, 24];     % Heater 1 (nodos 815 y 816)
idxH2 = 34;           % Heater 2

% Matriz B (2 columnas: una por cada controlador)
B = zeros(N, 2); 
B(idxH2, 2) = 1 / tmm.C(idxH2); 
B(idxH1, 1) = (1/2) * (1 ./ tmm.C(idxH1)); % Reparte la potencia del H1 entre los dos nodos

%matriz identidad para comprobar x
CI = eye(38);
DI = zeros(38,3);
DII = zeros(38,2);
% Matriz C_ss (2 filas: una para cada sensor/medida)
C_ss = zeros(2, N);                      
C_ss(1, idxH1) = 0.5;                     
C_ss(2, idxH2) = 1;

D_ss = zeros(2, 2);



T0 = tmm.T0;  % Temperatura setpoint inicial = T0

%  B original  para los heaters (u)
B_u = zeros(N, 2);
B_u(idxH1, 1) = 0.5 ./ tmm.C(idxH1);  % Repartir H1
B_u(idxH2, 2) = 1 ./ tmm.C(idxH2);     % H2 completo

% Creamos una B para las perturbaciones en el radiador(wd)
% Si wd entra en todos los nodos, Bd es la identidad dividida por la capacidad

idx_rad = 9:32;  % D801-D824 
B_d = zeros(N,1);
B_d(idx_rad,1) = (1/24) ./ tmm.C(idx_rad);


% Ahora la entrada U del bloque será un vector [u; wd]
B_total = [B_u, B_d]; 


% La matriz D debe  coincidir con las dos entradas
D_total = [D_ss, zeros(size(C_ss, 1), size(B_d, 2))];

sys = ss(A, B_total, C_ss, D_total);


%% 5. PI 
Kp1 = 60 ;    Ki1 = 0.5;
Kp2 = 2 ;    Ki2 = 0.5;

%% parámetros 
% ajustar según el caso a estudiar
Pmin = 5; Vmin = 28;
V0 = 30; 
R0 = (Vmin^2)/Pmin; 
P0 = (V0^2)/R0;
Vmax = 34;
Pmax = (Vmax^2)/R0;

dVdt = -4/60; % de 34 a 28 +-4v/min 
dVdt = 0;
Vb0 = 30;


TH2 = T0(34);
TH1 = (T0(23)+T0(24))/2;
Tsp = [TH1; TH2];

RH1 = R0;
RH2 = R0;

Q0_lin = [(Q0_linear(23)+Q0_linear(24)); Q0_linear(34)];
Ufeedforward = [(Q0(23)+Q0(24)) ; Q0(34)];

step_temp = 1; %o 2 según caso ensayado
steptime = 20;

wd = 0;



%% ruido
sigma = 0.0025;


%% sistema aumentado

ny = size(C_ss,1);

% continuo
Ts = 0;
A_aug = [A zeros(N,ny);
        -C_ss zeros(ny)];

B_aug = [B_u;
         zeros(ny,2)];

Ts = 0;

% discretización
% Ts = 1;
% sys_d = c2d(sys,   Ts, 'zoh');
% 
% [Ad, Bd, Cd, Dd] = ssdata(sys_d);  
% 
% sys_d_full = c2d(ss(A, [B_u, B_d], C_ss, zeros(ny,3)), Ts, 'zoh');
% Bd_u  = sys_d_full.B(:,1:2);   % parte de heaters
% Bd_d  = sys_d_full.B(:,3);     % parte de perturbación
% 
% 
% Ad_aug = [Ad,          zeros(N, ny);
%           -Cd*Ts,      eye(ny)     ];
% 
% Bd_aug = [Bd_u;
%           zeros(ny,2) ];



%% LQR

%continuo
% Q_lqr = diag([ones(1,N)*0.001, 1, 10]);
% R_lqr = diag([100, 50]);
Q_x = ones(1,N) * 0.001;
Q_x(idxH1) = 0.5;    
Q_x(idxH2) = 0.25;
Q_lqr = diag([Q_x, 20, 10]); 

% Q_lqr(idxH1, idxH1) = diag([1, 1]);          
% Q_lqr(idxH2, idxH2) = 1;                     

R_lqr = diag([3 , 12 ]);

% %discretizado 
% K_total = dlqr(Ad_aug, Bd_aug, Q_lqr, R_lqr);
% K    = K_total(:, 1:N);       % ganancia de estado
% Kint = K_total(:, N+1:N+ny);  % ganancia integral


%continuo

K_total = lqr(A_aug, B_aug, Q_lqr, R_lqr);
K = K_total(:,1:N);
Kint = K_total(:,N+1:N+2);

%% Kalman

Qnc = 1e-1*eye(N);
Qnc(idx_rad,    idx_rad)    = 6.25e-2*eye(24);
Qnc(idxH1(1),    idxH1(1))    = 1e-5;
Qnc(idxH1(2),    idxH1(2))    = 1e-5;
Qnc(idxH2,    idxH2)    = 1e-5;

Rnc = 6.25e-6 * eye(2);

Qn = Qnc;
Rn = Rnc;

% continuo
Kf = lqe(A, eye(N),C_ss, Qn, Rn );
syskf= ss(A, B_u, C_ss, D_ss);


% % discretizado
% 
% Qn = Qnc.*Ts;
% Rn = Rnc./Ts;
% 
% [P_kf_d, ~, Kf_d_raw] = dare(Ad', Cd', Qn, Rn);
% 
% Kf = (P_kf_d * Cd') / (Cd * P_kf_d * Cd' + Rn);
% syskf= ss(Ad, Bd_u, Cd, D_ss, Ts);








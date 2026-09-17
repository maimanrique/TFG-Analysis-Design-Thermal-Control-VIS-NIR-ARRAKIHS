%% 1. Read the ESATAN thermal model

% Modelo directamente leído de ESATAN (quitando los nodos inactivos que no
% se utilizan

esatan_model_red = readThermalModel_red("ARRAKIHS_ACASE2_reduced", "Modelo_reducido");

%% 2. Reduce (remove boundary nodes) and linearize

% Modelo de ESATAN procesado para crear KL y KR y quitar los nodos de
% contorno con temperatura fijada

tmm_red = processModel_red(esatan_model_red);



%% Los datos de los modelos provienen de un caso estacionario resuelto, así que tienen que cumplir la ecuacion de balance térmico (con una tolerancia)
%% Ecuación de balance completa
error1_red = tmm_red.KL*tmm_red.T0+tmm_red.KR*tmm_red.T0.^4+tmm_red.Q0;
max_error1_red = sum(abs(error1_red))

%% Ecuación de balance linealizada
error2_red = tmm_red.K0*tmm_red.T0+tmm_red.Q0_linear;
max_error2_red = sum(abs(error2_red))


%% 1. Read the ESATAN thermal model

% Modelo directamente leído de ESATAN (quitando los nodos inactivos que no
% se utilizan
esatan_model = readThermalModel("ARRAKIHS_ACASE2", "Thermal_models/Simple_model");


%% 2. Reduce (remove boundary nodes) and linearize

% Modelo de ESATAN procesado para crear KL y KR y quitar los nodos de
% contorno con temperatura fijada
tmm = processModel(esatan_model);






%% Los datos de los modelos provienen de un caso estacionario resuelto, así que tienen que cumplir la ecuacion de balance térmico (con una tolerancia)
%% Ecuación de balance completa
error1 = tmm.KL*tmm.T0+tmm.KR*tmm.T0.^4+tmm.Q0;
max_error1 = sum(abs(error1))

%% Ecuación de balance linealizada
error2 = tmm.K0*tmm.T0+tmm.Q0_linear;
max_error2 = sum(abs(error2))




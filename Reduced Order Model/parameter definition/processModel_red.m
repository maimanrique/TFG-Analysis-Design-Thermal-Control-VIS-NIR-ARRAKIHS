function out = processModel_red(tmm_red)
% processModel  Reduce and linearize a thermal model
%
%   out = processModel(tmm)
%
% Takes a struct produced by readThermalModel (with fields nD, nB, nN, GL,
% GR, T0, C, Q0) and returns a linearized, reduced-order model in which
% boundary nodes have been eliminated and radiative couplings have been
% linearized around the operating point T0.
%
% The resulting state equation (temperatures in Kelvin) is:
%
%   C .* dT/dt = K0 * T + Q0_linear
%
% Outputs (struct "out")
%   out.N         number of nodes (diffusion only, boundary removed)
%   out.KL         linear coupling matrix (nN x nN) [W/K]
%   out.KR         radiative coupling matrix (nN x nN) [W/K^4]
%   out.Q0         original heat loads summing boundary contribution  (nN x 1) [W]
%   out.C          capacitance vector   (nN x 1) [J/K]
%   out.K0         linearized system conductance matrix (nN x nN) [W/K]
%   out.Q0_linear  total linearized constant heat-load vector (nN x 1) [W]
%                  (includes Q0, and radiative linearization residuals)

nD_red = tmm_red.nD;
nN_red = tmm_red.nN;

iD_red = 1:nD_red;            % diffusion-node indices
iB_red = (nD_red+1):nN_red;       % boundary-node indices


% =====================================================================
% 1. Create full KL and KR matrices including boundary nodes
% =====================================================================
% GL and GR are full symmetric matrices (both halves populated, zero diagonal).
% The conductance-matrix diagonal is the negative sum of each full row.
KL_diag_red = -sum(tmm_red.GL, 2);             % (nN x 1)
KL_red = tmm_red.GL + diag(KL_diag_red);           % (nN x nN)

KR_diag_red = -sum(tmm_red.GR, 2);             % (nN x 1)
KR_red = tmm_red.GR + diag(KR_diag_red);           % (nN x nN)
stph_boltz = 5.670374419e-8;   % Stefan-Boltzmann constant [W/m^2/K^4]
KR_red = KR_red * stph_boltz;   % convert to W/K^4


% =====================================================================
% 2. Remove boundary nodes (shall be at the bottom)
% =====================================================================
QLb_red = KL_red(iD_red, iB_red) * tmm_red.T0(iB_red);   % boundary contribution to heat loads (nD x 1)
QRb_red = KR_red(iD_red, iB_red) * tmm_red.T0(iB_red).^4;   % boundary contribution to radiative heat loads (nD x 1)

KL_red = KL_red(iD_red, iD_red);   % (nD x nD)
KR_red = KR_red(iD_red, iD_red);   % (nD x nD)
T0_red = tmm_red.T0(iD_red);   % (nD x 1)
C_red = tmm_red.C(iD_red);    % (nD x 1)
Q0_red = tmm_red.Q0(iD_red) + QLb_red + QRb_red;   % include boundary contribution to heat loads (nD x 1)

% =====================================================================
% 3. Linearize radiative couplings around T0
% =====================================================================
% From equations (4.16) and (4.17):
%
%   K_ij = KL_ij + 4 * KR_ij * T0_j^3          (4.16)
%   Q_i  = Q0_i  - 3 * sum_j( KR_ij * T0_j^4 ) (4.17)
%
% So the linearized system is  C .* dT/dt = K0 * T + Q0_linear

K0_red       = KL_red + 4 * KR_red * diag(T0_red.^3);       % (nD x nD) linearized conductance
Q0_linear_red = Q0_red - 3 * KR_red * (T0_red.^4);           % (nD x 1)  linearized heat load

% =====================================================================
% Pack outputs
% =====================================================================
out.N         = nD_red;
out.KL        = KL_red;
out.KR        = KR_red;
out.C         = C_red;
out.T0        = T0_red;
out.K0        = K0_red;
out.Q0        = Q0_red;
out.Q0_linear = Q0_linear_red;

end

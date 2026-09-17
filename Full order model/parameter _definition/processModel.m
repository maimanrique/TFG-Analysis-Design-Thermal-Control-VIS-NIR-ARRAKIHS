function out = processModel(tmm)
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

nD = tmm.nD;
nN = tmm.nN;

iD = 1:nD;            % diffusion-node indices
iB = (nD+1):nN;       % boundary-node indices


% =====================================================================
% 1. Create full KL and KR matrices including boundary nodes
% =====================================================================
% GL and GR are full symmetric matrices (both halves populated, zero diagonal).
% The conductance-matrix diagonal is the negative sum of each full row.
KL_diag = -sum(tmm.GL, 2);             % (nN x 1)
KL = tmm.GL + diag(KL_diag);           % (nN x nN)

KR_diag = -sum(tmm.GR, 2);             % (nN x 1)
KR = tmm.GR + diag(KR_diag);           % (nN x nN)
stph_boltz = 5.670374419e-8;   % Stefan-Boltzmann constant [W/m^2/K^4]
KR = KR * stph_boltz;   % convert to W/K^4


% =====================================================================
% 2. Remove boundary nodes (shall be at the bottom)
% =====================================================================
QLb = KL(iD, iB) * tmm.T0(iB);   % boundary contribution to heat loads (nD x 1)
QRb = KR(iD, iB) * tmm.T0(iB).^4;   % boundary contribution to radiative heat loads (nD x 1)

KL = KL(iD, iD);   % (nD x nD)
KR = KR(iD, iD);   % (nD x nD)
T0 = tmm.T0(iD);   % (nD x 1)
C = tmm.C(iD);    % (nD x 1)
Q0 = tmm.Q0(iD) + QLb + QRb;   % include boundary contribution to heat loads (nD x 1)

% =====================================================================
% 3. Linearize radiative couplings around T0
% =====================================================================
% From equations (4.16) and (4.17):
%
%   K_ij = KL_ij + 4 * KR_ij * T0_j^3          (4.16)
%   Q_i  = Q0_i  - 3 * sum_j( KR_ij * T0_j^4 ) (4.17)
%
% So the linearized system is  C .* dT/dt = K0 * T + Q0_linear

K0       = KL + 4 * KR * diag(T0.^3);       % (nD x nD) linearized conductance
Q0_linear = Q0 - 3 * KR * (T0.^4);           % (nD x 1)  linearized heat load

% =====================================================================
% Pack outputs
% =====================================================================
out.N         = nD;
out.KL        = KL;
out.KR        = KR;
out.C         = C;
out.T0        = T0;
out.K0        = K0;
out.Q0        = Q0;
out.Q0_linear = Q0_linear;

end

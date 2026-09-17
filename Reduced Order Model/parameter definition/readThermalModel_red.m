function out = readThermalModel_red(model_name, folder)
% readThermalModel  Read ESATAN thermal model files (.nl.csv, .nd.csv, .gl.csv, .gr.csv)
%
%   out = readThermalModel(model_name)
%   out = readThermalModel(model_name, folder)
%
% Inputs
%   model_name : string/char, base name shared by all model files
%   folder     : optional, folder containing the files (default: pwd)
%
% Outputs (struct "out")
%   out.nD          number of D (diffusion) nodes
%   out.nB          number of B (boundary) nodes
%   out.nN          total active nodes (nD + nB)
%   out.GL          GLs matrix (nN x nN), linear conductances
%   out.GR          GR  matrix (nN x nN), radiative conductances
%   out.T0           temperature vector (nN x 1)
%   out.C           capacitance vector  (nN x 1)
%   out.Q0           heat load vector    (nN x 1)
%
% Notes
% - Ignores comment/header lines starting with '#'.
% - Node list lines (.nl.csv) are CSV: Type,Number,"Label",Model
% - Node data (.nd.csv) has columns: T, C, Q (one row per node).
% - GLs and GR files have header lines then numeric CSV rows.
% - Inactive (X) nodes are removed from all outputs.
% - All B nodes must appear after all D nodes; an error is raised otherwise.

if nargin < 2 || isempty(folder)
    folder = pwd;
end

nlFile = fullfile(folder, model_name+'.nl.csv');
ndFile = fullfile(folder, model_name+'.nd.csv');
glFile = fullfile(folder, model_name+'.gl.csv');
grFile = fullfile(folder, model_name+'.gr.csv');

% ---- Read node list ----
typeChar = readNodeTypes(nlFile);
Nfull    = numel(typeChar);

maskX = (typeChar == 'X');
keepIdx = find(~maskX);

% ---- Remove X nodes and reorder type vector ----
typeActive = typeChar(keepIdx);

% ---- Validate ordering: all B nodes must be at the end ----
maskD_active = (typeActive == 'D');
maskB_active = (typeActive == 'B');

nD = nnz(maskD_active);
nB = nnz(maskB_active);
nN = nD + nB;

if nN ~= numel(typeActive)
    error('Unknown node types found (only D and B expected after removing X nodes).');
end

lastD = find(maskD_active, 1, 'last');
firstB = find(maskB_active, 1, 'first');
if ~isempty(lastD) && ~isempty(firstB) && (firstB < lastD)
    error('Invalid node ordering: B nodes must all appear after D nodes.');
end

% ---- Read node data (T, C, Q) from .nd.csv ----
ndData = readNodeData(ndFile, Nfull);
ndData = ndData(keepIdx, :);

% ---- Read GLs matrix (GL) ----
GL_full = readConductanceMatrix(glFile, Nfull);
GL = GL_full(keepIdx, keepIdx);

% ---- Read GR matrix (GR) ----
GR_full = readConductanceMatrix(grFile, Nfull);
GR = GR_full(keepIdx, keepIdx);

% ---- Pack outputs ----
out.nD = nD;
out.nB = nB;
out.nN = nN;
out.GL = GL;
out.GR = GR;
out.T0  = ndData(:, 1) + 273.15;   % convert from Celsius to Kelvin
out.C  = ndData(:, 2);
out.Q0  = ndData(:, 3);

end

% ======================================================================

function typeChar = readNodeTypes(nlFile)
% Read node types from .nl.csv (D/B/X per row).
fid = fopen(nlFile, 'rt');
if fid < 0
    error('Cannot open node list file: %s', nlFile);
end
cleanup = onCleanup(@() fclose(fid));

types = char.empty(0,1);

while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end
    if startsWith(line, '#'), continue; end

    % First field is the node type letter
    % Example: D,500,"Strap1",ARRAKIHS_ACASE2
    k = find(line == ',', 1, 'first');
    if isempty(k)
        continue;
    end
    t = strtrim(line(1:k-1));
    if ~isempty(t)
        types(end+1,1) = t(1); %#ok<AGROW>
    end
end

typeChar = types;
end

% ======================================================================

function data = readNodeData(ndFile, Nfull)
% Read node data (T, C, Q) from .nd.csv.  Returns Nfull x 3 matrix.
fid = fopen(ndFile, 'rt');
if fid < 0
    error('Cannot open node data file: %s', ndFile);
end
cleanup = onCleanup(@() fclose(fid));

rows = cell(0,1);

while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end
    if startsWith(line, '#'), continue; end

    nums = parseNumericCSV(line);
    if isempty(nums), continue; end

    rows{end+1,1} = nums; %#ok<AGROW>
end

if numel(rows) < Nfull
    error('Node data file has only %d rows, expected at least %d.', numel(rows), Nfull);
end

% Take the last Nfull rows
data = zeros(Nfull, 3);
for i = 1:Nfull
    v = rows{end - Nfull + i};
    data(i, 1:min(3, numel(v))) = v(1:min(3, numel(v)));
end
end

% ======================================================================

function M = readConductanceMatrix(matFile, Nfull)
% Read a full numeric conductance matrix (GLs or GR) from CSV,
% skipping header/comment lines. Returns Nfull x Nfull.
fid = fopen(matFile, 'rt');
if fid < 0
    error('Cannot open conductance file: %s', matFile);
end
cleanup = onCleanup(@() fclose(fid));

rows = cell(0,1);
maxCols = 0;

while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end
    if startsWith(line, '#'), continue; end

    nums = parseNumericCSV(line);
    if isempty(nums), continue; end

    rows{end+1,1} = nums; %#ok<AGROW>
    maxCols = max(maxCols, numel(nums));
end

if isempty(rows)
    error('No numeric data found in file: %s', matFile);
end

nRows = numel(rows);
if nRows < Nfull
    error('File %s has only %d numeric rows, expected at least %d.', matFile, nRows, Nfull);
end

% Build numeric array and take last Nfull rows / columns
A = NaN(nRows, maxCols);
for i = 1:nRows
    v = rows{i};
    A(i, 1:numel(v)) = v;
end

A = A(end-Nfull+1:end, :);

if size(A,2) < Nfull
    error('File %s has only %d columns, expected at least %d.', matFile, size(A,2), Nfull);
end

M = A(:, end-Nfull+1:end);
M(isnan(M)) = 0;
end

% ======================================================================

function nums = parseNumericCSV(line)
% Parse a CSV line of numbers, tolerant to spaces and trailing commas.
% Returns double row vector (empty if not numeric).
line = regexprep(line, ',\s*$', ''); % remove trailing comma(s)
parts = regexp(line, '\s*,\s*', 'split');

nums = zeros(1, numel(parts));
for k = 1:numel(parts)
    x = str2double(parts{k});
    if isnan(x)
        nums = []; % not a numeric row
        return;
    end
    nums(k) = x;
end
end

function modeli = computeEms(networkFile,emComputationTool)
% required inputs
% - networkFile: network file name
%
% optional inputs
% - emComputationTool: software tool for em computation  
%    - 'metatool' (by default)
%    - 'efmtool'
%
% outputs
% - modeli
%   - .m: names of intracellular metabolites 
%   - .x: names of extracellular metabolites 
%   - .r: names of reactions
%   - .sm: stoichiometric matrix for m
%   - .sx: stoichiometric matrix for x
%   - .z: the full set of ems
%   - .sxz: sx*z

clc
% check inputs
if nargin < 2
    emComputationTool = 'metatool'; % default
end

% check if network file exists
icheck = exist(networkFile,'file');
if icheck == 0
    error('Network file ''%s'' does not exist!',networkFile)
end

modeli = runEmSoftwaretool(networkFile,emComputationTool);
save('myPremodels.mat','modeli')
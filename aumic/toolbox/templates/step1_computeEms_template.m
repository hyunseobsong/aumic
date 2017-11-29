clear all
%--------------------------------------------------------------------------
networkFile = ''; % network file name (e.g., 'network.txt')
emComputationTool = 'efmtool'; % 'efmtool' (default) or 'metatool'
%==========================================================================

modeli = computeEms(networkFile,emComputationTool)
edit step2_processEms

% networkFile
% - metabolic network file 
% - example: 
%   - networkFile = 'yeastNetwork.txt'; 

% emComputationTool
% - software tool for computing EMs from the network
% - 'efmtool' (default) or 'metabool'
% - example:
%   - emComputationTool = 'efmtool';
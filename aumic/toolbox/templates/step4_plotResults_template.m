%--------------------------------------------------------------------------
exportResults2excel = 'yes'; % 'yes' or 'no'
vars2plot = {
    'mets'
    }; % variables to plot ('mets','fluxes','rz','u','v','erel')
mets4plot = {
    }; % specific metabolites to plot
ems4plot = [
    ]; % specific EMs to plot 'rz', 'u', 'v', and/or 'erel'
fluxes4plot = {
    }; % speific fluxes to plot 
%==========================================================================

plotResults(result,model,vars2plot,mets4plot,fluxes4plot,...
    ems4plot,exportResults2excel)

% vars2plot
% - variables to plot
% - 'mets', 'fluxes', 'rz', 'u', 'v', 'erel', 'vol'
% - example:
%     vars2plot = {
%         'mets','fluxes','u','v'
%         };

% mets4plot
% - list of metabolites to plot
% - example:
%     mets4plot = {
%         'GLCx','XYL','vol'
%         };
% ems4plot
% - list of EMs to plot
% - example:
%     ems4plot = [
%         3,1
%         ];

% fluxes4plot
% - list of fluxes to plot
% - example:
%     fluxes4plot = {
%         'R1','R31'
%         };

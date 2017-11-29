%--------------------------------------------------------------------------
famRule = {
    }; % 1) EM classification rule, 2) reference flux for EM normalization
met.biom = ''; % metabolite name for biomass
emProcess.how = ''; % 'reduce' or 'lump' or 'none'
emProcess.data.fam{1} = {
    }; % exp data (relative fluxes) for EM reduction or lumping 
%==========================================================================

model = processEms(modeli,famRule,met,emProcess)
edit step3a_simulate
edit udf_basic
edit step3a_simulate

% famRule: 
% - rules for EM classification into families
% - each row defines each family 
% - first column: rule (make rules using 'and', 'or', 'not')
% - second column: reference flux for normalization 
% - example: 
%     famRule = { ...
%         'R1 and (not R31)','R1'
%         '(not R1) and R31','R31'
%         };

% met.biom: 
% - metabolite name for biomass
% - example: met.biom = 'BIOM';

% emProcess.how: method for EM processing ('none' or 'reduce' or 'lump')
% - example: % emProcess.how = 'reduce'; 

% emProcess.data.fam{i}: 
% - experimental data used for EM processing for family i
% - not needed if emProcess.how = 'none'
% - first column: flux name
% - second column: flux value
% - example: 
%     emProcess.data.fam{1} = {
%         'R1', -1 % GLC consumption
%         'R35', 0.0169 % BIOM production
%         'R12b', 1.8878 % ETHx production
%         'R6', 0.0556 % GOLx production
%         };
%     emProcess.data.fam{2} = {
%         'R31', -1 % XYL consumption 
%         'R35', 0.0249 % BIOM production 
%         'R12b', 1.1423 % ETHx production 
%         'R6', 0.0474 % GOLx production 
%         'R32', 0.049 % XOLx production 
%         };

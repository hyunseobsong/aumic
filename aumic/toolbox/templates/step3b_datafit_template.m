%--------------------------------------------------------------------------
datafit.file = ''; % excel datafile name (e.g., 'datafile.xlsx')
datafit.tab = {
    }; % excel sheet name
datafit.initial.met = {
    }; % initial metabolite concentrations for each experiment
datafit.initial.relEnz = [
    ]; % initial (relative) enzyme levels for each experiment
datafit.para = {
    }; % parameters to optimize
datafit.displayFrequency = 10; % display frequency 
closeFigs = 'yes'; % 'yes' or 'no'
%==========================================================================

kinetic = modelfit(model,kinetic,cybernetic,datafit,initial,closeFigs);
[result,model] = simulate(model,initial,tspan,kinetic,cybernetic,options4sim)
edit step4_plotResults
% [kinetic,result,model] = modelfit(model,kinetic,cybernetic,datafit,...
%     tspan,initial,closeFigs,options4sim);
% edit step4_plotResults

% datafit.file 
% - experimental data
% - .file: excel file name
% - .tab: tab name
% - .initial.met: intial metabolic concentrations for each exp. 
%    - first column: metabolite names
%    - second column: initial concentration for the first exp.
%    - third column: initial concentration for for the second exp., etc.
% - datafit.initial.relEnz: initial enzyme concetnration for each exp.
%    - first column: initial concentration for the first exp.
%    - second column: initial concentration for for the second exp., etc.
% - datafit.para: parameters to optimize through data fitting
% - datafit.displayFrequency: frequency of displaying parameter values
% - example:
%     datafit.file = 'yeastData.xlsx';
%     datafit.tab = {
%         'glucose'
%         'xylose'
%         };
%     datafit.initial.met = {
%         'BIOM',0.1, 0.1
%         'GLCx',636.7, 0
%         'XYL',0, 347.45
%         'ETH',0,4.04
%         'GOLx',0,1
%         'XOLx',0 0
%         };
%     datafit.initial.relEnz = [
%         0.0909,0.0909
%         0.0909,0.0909
%         0.0909,0.0909
%         0.513,0.513
%         0.513,0.513
%         0.513,0.513
%         0.513,0.513            
%         ];
%     datafit.para = {
%         'kmaxG1'
%         'kmaxG2'
%         'kmaxG3'
%         'kmaxX1'
%         'kmaxX2'
%         'kmaxX3'
%         'kmaxX4'
%         'KinhG'
%         'KinhX'
%         };
%     datafit.displayFrequency = 10;

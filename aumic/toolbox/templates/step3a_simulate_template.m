%--------------------------------------------------------------------------
kinetic.kmax = { 
    }; % maximum reaction rate constants 
kinetic.K = {
    }; % all other kinetic parameters than kmax
initial.met = {
    }; % initial concentration of extracellular metabolites

cybernetic.metabObj= 'carbon uptake'; % 'carbon uptake' or 'growth'
cybernetic.ncarbon = [
    ]; % the number of carbon elements in the nutrients 

cybernetic.al = [
    0.01
    ]; % constitutive rate of enzyme synthesis
cybernetic.be = [
    0.02
    ]; % rate of enzyme degradation
cybernetic.ke = [
    1
    ]; % rate of enzyme degradation
initial.relEnz = [
    ]; % initial relative enzyme level
tspan = [
    ]; % time span for the integration of ODEs

cybernetic.vars2fix = {
    }; % cybernetic variables to fix to one ('u','v','erel')

options4sim.when2switch = [
   ]; % times at which switching occurs
options4sim.met2switch = {
    }; % metabolic concentrations at switching times 

wantDatafit = 'no'; % 'yes' or 'no'
%==========================================================================

[result,model] = simulate(model,initial,tspan,kinetic,cybernetic,options4sim)
switch wantDatafit
    case 'yes'
        edit step3b_datafit 
    case 'no'
        edit step3b_datafit 
        edit step4_plotResults
    otherwise
        error('wantDatafit ''%s'' is not a supported option!',wantDatafit)
end

% kinetic.kmax
% - maximum reaction rate constants
% - each row is associated with each EM 
% - first column: name of kmax's (as defined by user)
% - second column: numerical values
% - example:
%     kinetic.kmax = { 
%         'kmaxG1', 173
%         'kmaxG2', 32.1
%         'kmaxG3', 28.6
%         'kmaxX1', 18.8
%         'kmaxX2', 21
%         'kmaxX3', 14.3
%         'kmaxX4', 11.9
%         }; 

% kinetic.K
% - all other kinetic parameters than kmax
% - each row is NOT necessarily associated with each EM 
% - first column: name of K's (as defined by user)
% - second column: numerical values
% - example
%     kinetic.K = { 
%         'KG1', 3.14
%         'KG2', 3.14
%         'KG3', 3.14
%         'KX1', 22.7
%         'KX2', 22.7
%         'KX3', 22.7
%         'KX4', 22.7
%         'KinhG', 410
%         'KinhX', 224
%         };  

% initial.met
% - initial concentration of (extracellular) metabolites
% - 'model.x' to get the list of metabolites
% - provide nonzero concentrations only (to save space)
% - first column: name of metabolites
% - second column: numerical values
% - example:
%     initial.met = {
%         'BIOM',2.3
%         'GLCx',444.44
%         'XYL',266.67
%         }; 

% cybernetic
% - cybernetic model-associated parameters
% - .metaObj: metabolic objective ('growth' or 'carbon uptake'), scalar
% - .ncarbon: number of carbons consumed through each EM, nz-by-1 
%            (not needed if metaObj = 'growth')
% - .al: constitutive enzyme synthesis rates, nz-by-1 or scalar
% - .be: enzyme degradation rates, nz-by-1 or scalar
% - .ke: reaction constant of inductive enzyme synthesis rates, nz-by-1 or scalar 
% - example:
%     cybernetic.metabObj= 'carbon uptake'; % 'carbon uptake' or 'growth'
%     cybernetic.ncarbon = [6 6 6 5 5 5 5]; 
%     cybernetic.al = 0.1;
%     cybernetic.be = 0.2;
%     cybernetic.ke = 1;

% initial.relEnz
% - initial concentrations of relative enzyme, nz-by-1
% - example:
%     initial.relEnz = [
%         0.0909
%         0.0909
%         0.0909
%         0.513
%         0.513
%         0.513
%         0.513
%         ]; 

% tspan
% - time span for integartion of ODEs
% - example: 
%     tspan = [0 40];

% cybernetic.fixedVars
% - cybernetic variabls to be forced to constant (i.e., 1)
% - effective for simulating non-cybernetic dynamic models for comparison
% - example:
%     cybernetic.fixedVars = {
%         'u'
%         'v'
%         'erel'
%         };

% wantDatafit: option for data fitting
% - choose either 'yes' or 'no'


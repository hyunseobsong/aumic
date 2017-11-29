% function [kinetic,result,model] = modelfit(model,kinetic,cybernetic,...
%     datafit,tspan,initial,closeFigs,options4sim)
function kinetic = modelfit(model,kinetic,cybernetic,...
    datafit,initial,closeFigs)

global iter4fit
iter4fit = 0;

if strcmp(closeFigs,'yes')
    close all
end

% get initial relative enzyme level from simulate.m setting
if isempty(datafit.initial.relEnz)
    datafit.initial.relEnz = repmat(initial.relEnz,1,length(datafit.tab));
end


% check datafit.file name
if exist(datafit.file,'file') == 0;
    error('datafit.file ''%s'' does not exist!',datafit.file)
end
datafit.tab = datafit.tab(:);

ntab = length(datafit.tab);

% check datafit.initial.met
if ~isempty(datafit.initial.met)
    ncol = size(datafit.initial.met,2);
    if ncol~=ntab+1
        error('The number of columns of datafit.initial.met should be %d.',ntab+1)
    end
end

% check datafit.initial.relEnz
nz = size(model.z,2);
if ~isempty(datafit.initial.relEnz)
%     error('Data4fit.initial.relEnz should not be empty!')
% else
    [nrow,ncol] = size(datafit.initial.relEnz);
    if nrow~=nz || ncol~=ntab
        error('The dimension of datafit.initial.relEnz should be %d-by-%d.',nz,ntab)
    end
end

% INITIAL METABOLITE CONCCENTRATIONS 
% initialize 
initialmet_ =  [model.x num2cell(zeros(length(model.x),ntab))];
% read initial condition from excel sheet
for itab = 1:ntab
    [data(itab).num,data(itab).txt] = xlsread(datafit.file,datafit.tab{itab});
    xname = data(itab).txt;
    xvalue = data(itab).num(1,2:end);
    xname(1) = [];
    iixdata = [];
    for ixdata = 1:length(xname)
        ixdata_ = find(strcmp(xname{ixdata},model.x));
        iixdata = [iixdata;ixdata_];
    end
    initialmet_(iixdata,1+itab) = num2cell(xvalue');
end

% update intial condition from user input (i.e., datafit.initial.met)
if ~isempty(datafit.initial.met)
	xname = datafit.initial.met(:,1);
    iixdata = [];
    for ixdata = 1:length(xname)
        ixdata_ = find(strcmp(xname{ixdata},model.x));
        iixdata = [iixdata;ixdata_];
    end    
    for itab = 1:ntab
        initialmet_(iixdata,1+itab) = datafit.initial.met(:,1+itab);
    end
end
datafit.initial.met = initialmet_;

% % check for other settings
% for itab = 1:ntab
%     initial.met = datafit.initial.met;
%     initial.relEnz = datafit.initial.relEnz(:,itab);
%     cybernetic = checkSimulationSetting(model,initial,kinetic,cybernetic);
% end
cybernetic = checkSimulationSetting(model,initial,kinetic,cybernetic);

kmax.num = cell2mat(kinetic.kmax(:,2));
kmax.name = kinetic.kmax(:,1);
K.num = cell2mat(kinetic.K(:,2));
K.name = kinetic.K(:,1);

% check para2fit
para = datafit.para;
para = para(:);

np = length(para);
for ip = 1:np
    ip_check = find(strcmp(para{ip},[kmax.name;K.name]));
    if isempty(ip_check)
        error('para2fit name ''%s'' is not found!',para{ip})
    end
end

% search kmax
iip.kmax = []; ip_ = [];
for ip = 1:np
    ip_ = find(strcmp(para{ip},kmax.name));
    iip.kmax = [iip.kmax;ip_];
end

% search K
iip.K= []; ip_check = [];
for ip = 1:np
    ip_ = find(strcmp(para{ip},K.name));
    iip.K = [iip.K;ip_];
end

% collect info for datafit
modeltemp = model;
modeltemp.datafit.initial = datafit.initial;
modeltemp.kinetic = kinetic;
modeltemp.cybernetic = cybernetic;
modeltemp.displayFrequency = datafit.displayFrequency;
if isempty(modeltemp.displayFrequency)
    modeltemp.displayFrequency = 10; % default value    
end

modeltemp.tab = datafit.tab;
modeltemp.ntab = ntab;
modeltemp.data = data;
modeltemp.para = para;
modeltemp.iip = iip;
modeltemp.fitstatus='notyet';
modeltemp.statistics='notyet';

% confirm user inputs
row1{1} = 'Met';
for itab = 1:ntab
    row1{1,1+itab} = ['exp',num2str(itab)];
end
confirmMet = [row1;datafit.initial.met];

row1{1} = 'Enzyme';
for iz=1:nz
    col1{iz,1} = ['e_rel,',num2str(iz)'];
end
confirmRelEnz = [row1;col1, num2cell(datafit.initial.relEnz)];
confirmMet
confirmRelEnz

% initial guesses
p0 = [kmax.num(iip.kmax);K.num(iip.K)];
popt = p0;
options = [];
[popt,fval,exitflag]=fminsearch(@paraoptObj,p0,options,modeltemp);
% % Statistical analysis
% model.statistics='doitnow';
% getStatistics(popt,model,fval,iter4fit)

% Draw smooth curves using paraoptObj
modeltemp.fitstatus='final';
paraoptObj(popt,modeltemp);

% optimal paramtere values:
disp(' ')
disp('Optimal parameter values:')
for ip=1:length(popt)
    fprintf('%10s %12.4e \n',char(modeltemp.para(ip)),popt(ip))
end

% get R2 values
modeltemp.fitstatus='notyet';
[f,R2,R2_adjusted] = paraoptObj(popt,modeltemp);
R2
R2_adjusted

kinetic = updateParameters(kinetic,datafit.para,popt);
% [result,model] = simulate(model,initial,tspan,kinetic,cybernetic,options4sim);
% edit step4_plotResults

function [zInFam,eta] = lumpEms(z_,iData,dataValue,model)

% global model

model.iData = iData;
model.dataValue = dataValue;
model.z_ = z_;

if isempty(iData)
    model.neta = 10; 
    a0 = [];
    eta = getEta(model,iData,a0);
    zInFam = z_*eta/sum(eta);
    sprintf('EM lumping without data... done.')
    return
end


options = [];
a0 =zeros(length(dataValue),1);
netamin = 1; netamax = 20; delneta = 1;
ineta = 0; AOPT = []; FVAL = []; NETA = [];
for neta = netamin:delneta:netamax
    ineta = ineta+1;
    model.neta = neta;
    
%     options = [];
%     options = gaoptimset(options,'Generations',100,'CrossoverFraction', 0.6);
%     rand('twister', 71); % These two commands are only included to
%     randn('state', 59); % make the results reproducible.
%     nvars = length(a0);
%     [a0, fval, exitflag, output, final_pop] = ga(@lumpObj,nvars,[],[],[],[],[],[],[],options);
% 
%     disp('----------------------------------------------------------')
%     disp('fminsearch starts!')
%     [aopt,fval,exitflag]=fminsearch(@lumpObj,a0,options);
    [aopt,fval,exitflag]=fminsearch(@lumpObj,a0,options,model);
%     [aopt,fval,exitflag]=fminunc(@lumpObj,a0,options,model);
    AOPT(:,ineta) = aopt;
    FVAL(ineta,1) = fval;
    NETA(ineta,1) = neta;
end
isol = find(FVAL==min(FVAL));
asol = AOPT(:,isol);
netasol = NETA(isol);

% final estimation
model.neta = netasol;
eta = getEta(model,iData,asol); 
zInFam = z_*eta/sum(eta);
% mean squred "relative" error
err=zInFam(iData)-dataValue;
error_fitting = norm(err./dataValue)^2/length(dataValue);

disp(' ')
disp('a-values')
disp(asol)
disp(sprintf('neta: %s',num2str(netasol)))

disp(' ')
disp('data vs model')
disp([dataValue zInFam(iData)])
disp('Mean squared "relative" error')
disp(error_fitting)
disp('-------------------------------------')

%--------------------------------------------------------------------------
% function f = lumpObj(a)
function f = lumpObj(a,model)
% global model

neta = model.neta;
z_ = model.z_;
iData = model.iData;
dataValue = model.dataValue;

eta = getEta(model,iData,a); 
zInFam = z_*eta/sum(eta);

err = zInFam(iData)-dataValue;
f = norm(err./dataValue)^2/length(dataValue);

%--------------------------------------------------------------------------
function eta = getEta(model,iData,a)

z_ = model.z_;
model.sxz = model.sx*z_;
iBiom = find(strcmp(model.biom,model.x));
ybiom = model.sxz(iBiom,:);
ybiom = ybiom(:);
neta = model.neta;

if isempty(iData)
    add = zeros(length(ybiom),1);
else
    add = a'*z_(iData,:);
    add = add(:);
end

eta = (ybiom + add).^neta;
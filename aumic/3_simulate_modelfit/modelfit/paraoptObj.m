function [f,R2,R2_adjusted,nV,nT,TDATA,YDATA,YMODEL] = paraoptObj(p,model)
global iter4fit h_fig

% prevent parameters from being negative 
for ip = 1:length(p)
    if p(ip) < eps
        p(ip) = eps;
        f = inf;
        return
    end
end

% initialization
f=[];R2=[];R2_adjusted=[];
nV=0;nT=[];TDATA=[];YDATA=[];YMODEL=[];

datafit.initial = model.datafit.initial;
kinetic = model.kinetic;
cybernetic = model.cybernetic;
displayFrequency = model.displayFrequency;

ntab = model.ntab;
data = model.data;
iip = model.iip;
nx = length(model.x);

nkmax4fit = length(iip.kmax);
for ip = 1:nkmax4fit
    kinetic.kmax{iip.kmax(ip),2} = p(ip);
end

nK4fit = length(iip.K);
for ip = 1:nK4fit
    kinetic.K{iip.K(ip),2} = p(nkmax4fit+ip);
end

YDATA = [];
YMODEL = [];
nT = [];
nV = 0;

for itab = 1:ntab
    xdata = data(itab).txt;
    xdata(1) = []; % delete time
    nxdata = length(xdata);
    iixdata = []; ixdata = [];
    for ixdata = 1:nxdata
        ixdata_ = find(strcmp(xdata{ixdata},model.x));
        iixdata = [iixdata;ixdata_];
    end
    
    initial.met = [datafit.initial.met(:,1),datafit.initial.met(:,1+itab)];
    initial.relEnz = datafit.initial.relEnz(:,itab);
    
    switch model.fitstatus
        case 'notyet'
            tspan = data(itab).num(:,1);
        case 'final'
            tspan = [data(itab).num(1,1),data(itab).num(end,1)]';
    end

    met.num = cell2mat(initial.met(:,2));
    met.name = initial.met(:,1);
    nx = length(model.x);
    x0 = zeros(nx,1);
    imet_ = [];
    for imet = 1:length(met.num)
        imet_ = find(strcmp(met.name{imet},model.x));
        x0(imet_) = met.num(imet);
    end

    erel0 = initial.relEnz(:);
    vol0 = 1;
    y0 = [x0;erel0;vol0];
    
    [result,model] = solveOdeModel(model,y0,tspan,kinetic,cybernetic);
    Tmodel = result.T;
    Ymodel = result.Y(:,iixdata);
    Tdata = data(itab).num(:,1);
    Ydata = data(itab).num(:,2:end);    
    result.DATA = Ydata;
    
    if strcmp(model.fitstatus,'notyet')
        % concatenate 
        [nt,nv]=size(Ydata);
        for j=1:nv
            nT = [nT;nt];
            for i=1:nt
                TDATA(i,nV+j)=Tdata(i);
                YDATA(i,nV+j)=Ydata(i,j);
                YMODEL(i,nV+j)=Ymodel(i,j);
            end
        end
        nV=nV+nv; % number of dependent variables
    end
    
%     if strcmp(model.statistics,'doitnow')
    if strcmp(model.fitstatus,'notyet')
        [nt,nv]=size(Ydata);
        for j=1:nv
            % Coefficient of determination (also referred to as R^2)
            SSerr = sum((Ydata(:,j)-Ymodel(:,j)).^2); %sum of squared residuals
            SSreg = sum((Ymodel(:,j)-mean(Ydata(:,j))).^2); %sum of squared residuals
            SStot = sum((Ydata(:,j)-mean(Ydata(:,j))).^2); 
            R2(j,itab) = 1 - SSerr/SStot; % R^2
            R2_adjusted(j,itab) = 1-(1-R2(j,itab))*(nt-1)/...
                (nt-length(p)-1); %Adjusted R^2
        end
    end
    
    if displayFrequency ~= inf
        if mod(iter4fit,displayFrequency) == 0
            if iter4fit == 0 
                h_fig(itab)=figure('name',model.tab{itab});
            else
                figure(h_fig(itab));
            end
            plotDatafit(Tmodel,result.Y(:,1:nx),model.x,Tdata,Ydata,iixdata)
        end
    end
    
end

if strcmp(model.fitstatus,'final'), return; end

if iter4fit == 0
    disp('Press any key to start data fitting.')
    pause
end

% YDATA
% nT
% nV
% pause

% Fitting error
% ---------------
% k=n_para_opt;
% v=nn_var;
% n=nn_YYexp;
% res=YYexp-YYmodel; % Residuals

% npara = length(p); % number of parameters
% nvar = nV; % number of variables 

k = length(p); % number of parameters
v = nV; % number of variables
n = nT; % a vector of the number of time points for each var
res=YDATA-YMODEL; % residual matrix

for m = 1:v
    % Sum of squared residuals for each dependent variable
    ssr(m) = res(1:n(m),m)'*res(1:n(m),m);
end

ssr=ssr(:); % ensuring column vector
s2 = ssr./(n-k/v); % Variances

% weighting factors
w_WSSRs = sum(n)./s2 / sum(n./s2); % Weights
for m=1:v
    w_WSSRm(m)=1/mean(YDATA(1:n(m),m))^2/n(m);
end
% w_WSSRm=w_WSSRm/sum(w_WSSRm)*v;
w_WSSRm=w_WSSRm/sum(w_WSSRm);
w_WSSRm=w_WSSRm(:);
w_UWSSR=ones(v,1);

% Measures for goodness of fit    
WSSRm = sum(w_WSSRm.*ssr); % Weighted sum of squared residuals (modified)
WSSRs = sum(w_WSSRs.*ssr); % Weighted sum of squared residuals (standard)
UWSSR = sum(w_UWSSR.*ssr); % Unweighted sum of squared residuals

f = WSSRm;

% if strcmp(paraopt_objfcn,'WSSRm')
%     f=WSSRm;
% elseif strcmp(paraopt_objfcn,'WSSRs')
%     f=WSSRs;
% elseif strcmp(paraopt_objfcn,'UWSSR')
%     f=UWSSR;
% end

% if initialchk4optPARA==1
%     initialchk4optPARA=0;
%     disp(' ')
%     display('Configure figures within 5 seconds!')
%     pause(5)
% end

if displayFrequency ~= inf
    if mod(iter4fit,displayFrequency)==0
        disp(' ')
        for j=1:k
            fprintf('%10s %12.4e \n',char(model.para(j)),p(j))
        end
        disp(' -----------------------------------------------------------------------')
        % disp(num2str(f))
        fprintf(' Iteration: %d\n',iter4fit)
        fprintf('Error(WSS): %10.4e\n',WSSRm)
    end
end
iter4fit = iter4fit + 1;

% % Display interim results
% disp(' ')
% if strcmp(paraopt_objfcn,'WSSRm')
%     fprintf('%6s %4d,%10s %11.4e,%6s %11.4e,%6s %11.4e\n',...
%         'Iter.',iter_paraOPT,'WSSRm(=J)',WSSRm,'WSSRs',WSSRs,'UWSSR',UWSSR)
% elseif strcmp(paraopt_objfcn,'WSSRs')
%     fprintf('%6s %4d,%6s %12.4e,%10s %11.4e,%6s %11.4e\n',...
%         'Iter.',iter_paraOPT,'WSSRm',WSSRm,'WSSRs(=J)',WSSRs,'UWSSR',UWSSR)
% elseif strcmp(paraopt_objfcn,'UWSSR')
%     fprintf('%6s %4d,%6s %12.4e,%6s %12.4e,%10s %11.4e\n',...
%         'Iter.',iter_paraOPT,'WSSRm',WSSRm,'WSSRs',WSSRs,'UWSSR(=J)',UWSSR)
% end
% disp(' -----------------------------------------------------------------------')
% for j=1:n_para_opt
%     fprintf('%10s %12.4e \n',char(para_opt_name(j)),para_opt(j))
% end
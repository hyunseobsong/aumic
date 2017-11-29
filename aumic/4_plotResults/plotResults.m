function plotResults(result,model,vars2plot,mets4plot,fluxes4plot,...
    ems4plot,exportResults2excel)

for ivar = 1:length(vars2plot)
    switch vars2plot{ivar}
        case 'mets'
            plotVars('mets',result,model,mets4plot)
        case 'fluxes'
            plotVars('fluxes',result,model,fluxes4plot)
        case 'rz'
            plotVars('rz',result,model,ems4plot)
        case 'u'
            plotVars('u',result,model,ems4plot)
        case 'v'
            plotVars('v',result,model,ems4plot)
        case 'erel'
            plotVars('erel',result,model,ems4plot)
        otherwise
            error('The variable ''%s'' to plot is not found!',vars2plot{ivar})
    end
end

switch exportResults2excel
    case 'yes'
        outputfile = 'results.xlsx';
        ifile = 0;
        while logical(exist(outputfile,'file'))
            ifile = ifile + 1;
            outputfile = ['results(',num2str(ifile),').xlsx'];
        end
        nx = length(model.x); 
        ny = size(result.Y,2);
        nz = size(model.z,2);
        
        TYhead = ['Time';model.x;'Vol']';
        TY = [result.T,result.Y(:,[1:nx,ny])];
        xlswrite(outputfile,TYhead,1,'A1')        
        xlswrite(outputfile,TY,1,'A2')      

        for iz=1:nz
            rzhead{iz,1} = ['rz',num2str(iz)];
            uhead{iz,1} = ['u',num2str(iz)];
            vhead{iz,1} = ['v',num2str(iz)];
            erelhead{iz,1} = ['erel',num2str(iz)];
            zhead{iz,1} = ['z',num2str(iz)];
        end
        TChead = ['Time';rzhead;uhead;vhead;erelhead]';
        TC = [result.T,result.RZ,result.U,result.V,result.EREL];
        xlswrite(outputfile,TChead,2,'A1')        
        xlswrite(outputfile,TC,2,'A2')   

        TFhead = ['Time';model.r]';
        TF = [result.T,result.FLUX];
        xlswrite(outputfile,TFhead,3,'A1')        
        xlswrite(outputfile,TF,3,'A2')   
        
        FZhead = ['Flux';zhead]';
        xlswrite(outputfile,FZhead,4,'A1')        
        xlswrite(outputfile,model.r,4,'A2')   
        xlswrite(outputfile,model.z,4,'B2')   
    case 'no'
    otherwise
        error('saveResult should be either ''yes'' or ''no''!')
        
end

function plotVars(plotwhat,result,model,uchoice,nmax)
% check inputs
if nargin < 5
    nmax = []; % the maximum number of subfigures 
    if nargin < 4
        uchoice = [];
    end
end

if isempty(nmax), nmax = 24; end

T = result.T;

switch plotwhat
    case 'mets'
        Y = result.Y;
        motherset = model.x;
        ny = size(result.Y,2);
        motherset{ny} = 'Vol';
        for ivar = 1:length(uchoice)
            ivar_ = find(strcmp(uchoice{ivar},motherset));
            if isempty(ivar_)
                error('The metabolite ''%s'' for plot is not found!',uchoice{ivar})
            end
        end
    case {'rz','u','v','erel'}
        nz = size(model.z,2);
        motherset_ = 1:nz;
        iiuchoice = [];
        for ivar = 1:length(uchoice)
            ivar_ = find(motherset_ == uchoice(ivar));
            if isempty(ivar_)
                error('The EM ''%d'' for plot does not exist!',uchoice(ivar))
            end
            iiuchoice = [iiuchoice;ivar_];
        end

        if strcmp(plotwhat,'rz')
            Y = result.RZ;
            for iz=1:nz
                motherset{iz,1} = ['r_{z,',num2str(iz),'}'];
            end
        elseif strcmp(plotwhat,'u')
            Y = result.U;
            for iz=1:nz
                motherset{iz,1} = ['u_{',num2str(iz),'}'];
            end
        elseif strcmp(plotwhat,'v')
            Y = result.V;
            for iz=1:nz
                motherset{iz,1} = ['v_{',num2str(iz),'}'];
            end
        elseif strcmp(plotwhat,'erel')
            Y = result.EREL;
            for iz=1:nz
                motherset{iz,1} = ['e_{rel,',num2str(iz),'}'];
            end
        end
        
        if ~isempty(iiuchoice)
            uchoice = motherset(iiuchoice);
        end
    case 'fluxes'
        Y = result.FLUX;
        motherset = model.r;
        for ivar = 1:length(uchoice)
            ivar_ = find(strcmp(uchoice{ivar},motherset));
            if isempty(ivar_)
                error('The flux ''%s'' for plot is not found!',uchoice{ivar})
            end
        end
%     case 'vol'
%         Y = result.Y(:,end);
%         motherset = {'vol'};
end

if isempty(uchoice),uchoice = motherset; end

allocateFigs(T,Y,motherset,uchoice,nmax)

function allocateFigs(T,Y,motherset,uchoice,nmax)

% indices of user-chosen vars
iivar = [];
for i = 1:length(uchoice)
    ivar = find(strcmp(uchoice{i},motherset));
    iivar = [iivar;ivar];
end

% check inputs
nvar = length(iivar);

nsub = ceil(nvar/nmax + eps);

for isub = 1:nsub
    if isub == nsub
        iirsub = iivar((isub-1)*nmax+1:end);
    else
        iirsub = iivar((isub-1)*nmax+1:isub*nmax);
    end
    varname = motherset(iirsub);
	plotSubfigs(T,Y(:,iirsub),iirsub,varname)
end

%--------------------------------------------------------------------------
function plotSubfigs(T,Y,iivar,varname)

nvar = length(iivar);

switch nvar
    case 1,nrow=1;ncol=1;
    case 2,nrow=2;ncol=1;
    case 3,nrow=3;ncol=1;
    case 4,nrow=2;ncol=2;
    case 5,nrow=3;ncol=2;
    case 6,nrow=3;ncol=2;
    case 7,nrow=4;ncol=2;
    case 8,nrow=4;ncol=2;
    case 9,nrow=5;ncol=2;
    case 10,nrow=5;ncol=2;
    case 11,nrow=4;ncol=3;
    case 12,nrow=4;ncol=3;
    case 13,nrow=5;ncol=3;
    case 14,nrow=5;ncol=3;
    case 15,nrow=5;ncol=3;
    case 16,nrow=4;ncol=4;
    case 17,nrow=5;ncol=4;
    case 18,nrow=5;ncol=4;
    case 19,nrow=5;ncol=4;
    case 20,nrow=5;ncol=4;
    case 21,nrow=6;ncol=4;
    case 22,nrow=6;ncol=4;
    case 23,nrow=6;ncol=4;
    case 24,nrow=6;ncol=4;
end

%--------------------------------------------------------------------------
figure()
for ivar=1:nvar
    subplot(nrow,ncol,ivar)
    plot(T,Y(:,ivar),'-');    
    ylabel(varname(ivar))    
    if nrow>=5
        set(gca,'xtick',[])
    end

    if ivar==nvar || ...
        (ncol==2 && ivar==nvar-1 && mod(nvar,ncol)==0) || ...
        (ncol==3 && ivar==nvar-1 && mod(nvar,ncol)==0) || ...
        (ncol==3 && ivar==nvar-2 && mod(nvar,ncol)==0) || ...
        (ncol==3 && ivar==nvar-1 && mod(nvar,ncol)==2) || ...
        (ncol==4 && ivar==nvar-1 && mod(nvar,ncol)==0) || ...
        (ncol==4 && ivar==nvar-2 && mod(nvar,ncol)==0) || ...
        (ncol==4 && ivar==nvar-3 && mod(nvar,ncol)==0) || ...
        (ncol==4 && ivar==nvar-1 && mod(nvar,ncol)==3) || ...
        (ncol==4 && ivar==nvar-2 && mod(nvar,ncol)==3) || ...
        (ncol==4 && ivar==nvar-1 && mod(nvar,ncol)==2)
        xlabel('Time')
        if nrow>=5
            set(gca,'xtickMode', 'auto')
        end
    end
end
drawnow

%--------------------------------------------------------------------------
function plotDatafit(T1,Y1,varname,T2,Y2,iixdata)

nvar = size(Y1,2);

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
for ivar=1:nvar
    subplot(nrow,ncol,ivar)
    plot(T1,Y1(:,ivar),'-');
    idata = find(iixdata==ivar);
    if ~isempty(idata)
        hold on
        plot(T2,Y2(:,idata),'o');
        hold off
        idata = idata + 1;
    end
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
	drawnow

end
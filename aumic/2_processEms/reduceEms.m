function zInFam = reduceEms(z,iData,dataValue,options)
if isempty(options)
    volRel = 0.99;
else
    volRel = options.volRel;
end

% sort (why? it sometimes helps to find solutions)
% -----------------------------------------------
zmas_row = z';
[zmas_row,iiSort]=sortrows(zmas_row,iData);
zmas = zmas_row';

ymas = zmas(iData,:);
nYmas = size(ymas,2);

%==========================================================================
% PHASE I: M_vertex
%==========================================================================
disp(' ')
disp('identification of vertices ...')
[kver,yver,vver] = getZver(ymas);
nYver = size(yver,1);
disp('... done')

%==========================================================================
% PHASE II: M_xx (e.g., M_99 or M_95)
%==========================================================================
if volRel<1.0
    disp(' ')
    disp('choice of main vertices ...')
    [KK_xx,K_xx,nEM_xx,V_xx_relative]=...
        get_M_xx(ymas',kver,vver,volRel);
%     Y_xx=ymas(K_xx,:);
    Y_xx=ymas(:,K_xx);
    disp('... done')
end

%==========================================================================
% PHASE III: M_act 
%==========================================================================
disp(' ')
disp('identification of active modes ...')

[mdataValue ndataValue]=size(dataValue);
if volRel < 1
    [K_active Y_active Weight error_fitting]=...
        get_M_act(K_xx,Y_xx',dataValue,ymas');
%     disp('   EM_xx     Yield Vector')
%     disp([K_xx' Y_xx])
else
    [K_active Y_active Weight error_fitting]=...
        get_M_act(kver,yver,dataValue,ymas');
end
Weight = Weight(:); % column vector
Y_active = Y_active';  
[mY_active]=size(Y_active,2);
disp('... done')

%--------
% Summary
%--------
if volRel==1.0
    disp(' ')
    disp(sprintf('    Mmas:     %d',nYmas))
    disp(sprintf('    Mver:     %d',nYver))
    disp(sprintf('    Mact:     %d',mY_active))
elseif volRel<1.0
    disp(' ')
    disp(sprintf('    Mmas:     %d',nYmas))
    disp(sprintf('    Mver:     %d',nYver))
    disp(sprintf('    Mxx:      %d',nEM_xx(end)))
    disp(sprintf('    Mact:     %d',mY_active))
    disp(sprintf('    Vxx/Vmas: %d',V_xx_relative(end)))
end

disp(' ')
disp('Active modes')
disp(Y_active)
disp('Weights')
disp(Weight)

Yestimated = Y_active*Weight;
disp('data vs model')
disp([dataValue Yestimated])
disp('Mean squared "relative" error')
disp(error_fitting)
disp('-------------------------------------')
zInFam = z(:,iiSort(K_active));

%==========================================================================
function [kver,yver,vver]=getZver(ymas)
ymas = ymas'; % row vector
nYmas = size(ymas,1);

% 'QJ' added by Philipp in June 2011
[KK_master V_master]=convhulln(ymas,{'Qt','QJ','Pp'}); 
[mKK_master nKK_master]=size(KK_master); 

%--------------------------------------------------------------------------
%Pick up yver
%--------------------------------------------------------------------------
ione=0;
for ii=1:nYmas 
    ij=0; 
    for i=1:mKK_master
        for j=1:nKK_master
            jj=KK_master(i,j);
            if ii==jj, ij=1; end
        end
    end
    if ij==1; ione=ione+1; kver(1,ione)=ii; end
end

yver=ymas(kver,:); 
[nYver nyver]=size(yver); mkver=length(kver);

%Compare V_master and vver 
[Kkver vver]=convhulln(yver,{'Qt','QJ','Pp'}); 
Vr=vver/V_master;
if abs(1-Vr)>=1.e-4
    disp('Note that vver/V_master is not unity!')
    disp('-------------------------------------')
    disp(Vr)
    pause
end

%==========================================================================
function [KK_xx,Kpenultimate,nEM_xx,V_xx_relative]=...
    get_M_xx(ymas,kver,vver,volRel)

[nYmas nymas]=size(ymas); 
yver=ymas(kver,:); 
[nYver nyver]=size(yver); mkver=length(kver);

% Find the mode maximizing each components
Ymax=zeros(nymas,1); iYmax=zeros(nymas,1);
for i=1:nYver 
    Yall=ymas(kver(i),:);
    for j=1:nyver
        if abs(Yall(j))>Ymax(j);
            iYmax(j)=kver(i);
            Ymax(j)=Yall(j);
        end
    end
end
mYmax=length(Ymax);

% disp('   EM         Ymax')
% disp('   ---------------')
% disp([iYmax Ymax])

%--------------------------------------------------------------------------
% Check Redundancy
%--------------------------------------------------------------------------
iYmaxsort=sort(iYmax); 
clear iYmax
ij=1; 
iYmax(ij,1)=iYmaxsort(1);
for i=2:mYmax 
    if iYmaxsort(i)~=iYmaxsort(i-1) 
        ij=ij+1;
        iYmax(ij,1)=iYmaxsort(i);
    end
end

%--------------------------------------------------------------------------
% Pick up the initial base EMs
%--------------------------------------------------------------------------

% Candidates for the initial base EMs
%-------------------------------------
Kmax=iYmax';%index for the modes maximizing BIOM and MAINT 
mKmax=length(Kmax);

%Knonmax: index for the modes other than Kmax
%----------------------------------------------
k=0;
for i=1:mkver
    iK=1;
    for j=1:mKmax
        if kver(i)==Kmax(j)
            iK=0;
        end
    end
    if iK==1
        k=k+1;
        Knonmax(1,k)=kver(i);
    end
end
mKnonmax=length(Knonmax); 

% Select Kbase using the largest distance 
%-------------------------------------
Kbase=Kmax; Knonbase=Knonmax;
mKbase=length(Kbase); mKnonbase=length(Knonmax);

for k=1:(nymas-mKmax+1)
    clear V p
    V=ymas(Kbase,:)';
    Distmax=eps;iDistmax=0;
    for i=1:mKnonbase
        p=ymas(Knonbase(i),:)';
        Dist=get_distance(V,p);
        if Dist>Distmax
            iDistmax=i; Distmax=Dist;
        end 
    end
    Kbase=[Kbase Knonbase(iDistmax)];
    mKbase=length(Kbase);
    Knonbase(iDistmax)=[];
    mKnonbase=length(Knonbase);

%     disp([mKbase Distmax])

end
format short

% disp('Compare Ymax and Ybase')
% disp('-----------------------')
% disp(Kmax)
% disp(Kbase)

[KKbase Vbase]=convhulln(ymas(Kbase,:),{'Qt','QJ','Pp'}); %Note 'Qt'
Vrbase=Vbase/vver;

% Expand EMs among the candidates 
%--------------------------------------
Kred=Kbase; Knonred=Knonbase;
mKred=length(Kbase); mKnonred=length(Knonred);

KK_xx=zeros(1,nYver);

% disp('    No of EMs V_xx/V_master')
% disp('    -----------------------------')

k=1; Vrmax=Vrbase;
KK_xx(k,1:mKbase)=Kred; V_xx_relative(k,1)=Vrmax; nEM_xx(k,1)=mKred;
%--------------------------------------------
%for k=2:mKnonbase+1
%--------------------------------------------
while Vrmax<=volRel; k=k+1;
%--------------------------------------------
    C=nchoosek(mKnonred,1); addset=nchoosek(Knonred,1); %addset=combnk(Knonred,1); 
    Vrmax=eps; iVrmax=0; 
    for i=1:C 
        Ksearch=[Kred,addset(i,:)];
        [KKsearch Vsearch]=convhulln(ymas(Ksearch,:),{'Qt','QJ','Pp'}); %Note 'Qt'
        Vr=Vsearch/vver;
        if Vr>Vrmax;
            iVrmax=i; Vrmax=Vr;
        end
    end
    Kred=[Kred, addset(iVrmax)];
    Knonred(iVrmax)=[];
    mKred=length(Kred); mKnonred=length(Knonred);
    KK_xx(k,1:mKred)=Kred;
    V_xx_relative(k,1)=Vrmax;
    nEM_xx(k,1)=mKred;
    
%     disp([mKred Vrmax])
    
end

Kpenultimate=sort(Kred);


%==========================================================================
function [K_active Y_active Weight error_fitting]=...
    get_M_act(K_xx,Y_xx,yieldData,ymas)

%Find nonzero component of yieldData: why? C and d are normalized w.r.t. yieldData
myieldData=length(yieldData); k=0; 
for i=1:myieldData 
    if yieldData(i)~=0 
        k=k+1; 
        iExpComp(k)=i; 
    end
end

% test ...
mK_xx=length(K_xx);

clear Kset
nKset=nchoosek(mK_xx,mK_xx); Kset=nchoosek(K_xx,mK_xx);

normmin=inf;
clear Ksetmin Cmin xmin
%=============================
for jj=1:nKset
%=============================
    %Define the input data to lsqlin
    clear C d C0 d0
    C=ymas(Kset(jj,:),iExpComp)'; d=yieldData(iExpComp)'; 
    %--------
    d=d(:);
    %--------
    C0=C; d0=d; [mC nC]=size(C);
    dtemp=repmat(d0,1,nC); C=C0./dtemp; d=ones(mC,1); %normalize C and d 
    Aeq=ones(1,nC); beq=1; %summation of x is unity.
    A=[];b=[]; %no inequality constraints
    lb=zeros(nC,1); ub=ones(nC,1); %lower and upper bound
    clear x0 x
    %x0=rand(nC,1);x0=x0/sum(x0);
    x0=repmat(1/nC,nC,1); %initial condition for x 
%     options=optimset('Display','off','MaxIter',100000); 
    options=optimset('Display','off'); 
%     [x,resnorm,residual,exitflag]=lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options);
    [x,resnorm,residual,exitflag]=lsqlin(C,d,A,b,Aeq,beq,lb,ub);
    
    x=abs(x); %eliminate components with minus values 

    if resnorm<normmin
        %inormmin=jj;
        Ksetmin=Kset(jj,:);
        normmin=resnorm;
        Cmin=C0;
        xmin=x;
    end

    %check if the solution is correct.
    if exitflag~=1
%         msgboxText{1} =  'lsqlin fails to find solution.'; 
%         msgbox(msgboxText,'lsqlin failure', 'warn');
        warning('The lsqlin solution might not be optimal.'); 
    end

%=============================
end
%=============================

%compare the prediction with exp data
Y1=Cmin*xmin;
Y2=d0;

k=0;
for i=1:nC
    if xmin(i)>1.e-6
        k=k+1;
        %K_active(k)=K_xx(i);
        K_active(k)=Ksetmin(i);
        Weight(k,1)=xmin(i);
    end
end

Y_active=ymas(K_active,:);

err=Y_active'*Weight-yieldData;

% mean squred "relative" error
error_fitting=norm(err./yieldData)^2/length(yieldData);

%==========================================================================
function dist=get_distance(V,p)

v0=V(:,end); v=V(:,1:end-1); [mv nv]=size(v);

for k=1:nv
    for i=1:nv
        A(k,i)=0;
        for j=1:mv
            A(k,i)=A(k,i)+(v(j,k)-v0(j))*(v(j,i)-v0(j));
        end
    end
end

for k=1:nv
    b(k,1)=0;        
    for j=1:mv
        b(k,1)=b(k,1)+(v(j,k)-v0(j))*(p(j)-v0(j));
    end
end

x=inv(A)*b;

sum=zeros(mv,1);
for i=1:nv
    sum=sum+x(i)*(v(:,i)-v0);
end
dist=norm(v0-p+sum);

%==========================================================================
function [G_row,y_BM,y_ATP_MNT,y_ATP_ANA] = ...
    get_GAR_equation(GAR,ir_anabolic,Sm_ATP,sxz,Z,ix_BM,ix_ATP)

% ATP consumption through EMs (for BM-producing and ATP-producing EM
% family)

y_BM=sxz(ix_BM,:);
y_ATP_MNT=sxz(ix_ATP,:);
y_ATP_ANA=(Sm_ATP(ir_anabolic)*Z(ir_anabolic,:));

if min(y_ATP_ANA) > 0
    error('ATP should not be produced through anabolic reactions!')
else
    y_ATP_ANA=abs(y_ATP_ANA);
end

G_row = y_ATP_MNT + y_ATP_ANA - GAR*y_BM;

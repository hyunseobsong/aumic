function dydt=odeModel(t,y,kin,cyb,model)

nx = length(model.x); ny = length(y); nz = size(model.z,2);
kmax = kin.kmax.num; al=cyb.al;be=cyb.be; ke=cyb.ke;
iBiom = strmatch(model.biom,model.x,'exact');

al=al(:); be=be(:); ke=ke(:);

dydt=zeros(ny,1);
for i=1:ny
    if y(i)<0,y(i)=eps;end
end
x=y(1:nx);
erel=y(nx+1:nx+nz);
vol = y(ny);

% rkin from a user-defined function (basic)
[rkin,flowin,flowout,xfeed,kla,xsat,dxdt_more] = udf_basic(t,kin,x,model);

rekin = zeros(nz,1);
for iZ=1:nz
    if kmax(iZ)==0
        rekin(iZ,1)=0;
    else
        rekin(iZ,1)=rkin(iZ)/kmax(iZ)*ke(iZ);
    end
end
% Return-On-Investment 
roi=erel.*rkin;
switch cyb.metabObj
    case 'carbon uptake'
        ncarbon=cyb.ncarbon;
        ncarbon=ncarbon(:);
        roi=ncarbon.*roi;
    case 'growth'
        roi=model.sxz(iBiom,:)'.*roi;
end

roi=max(roi,zeros(nz,1));
pu=max(roi,zeros(nz,1)); pv=max(roi,zeros(nz,1));
sumpu=sum(pu); maxpv=max(pv);
if sumpu>0, u=pu/sumpu; else u=zeros(nz,1); end
if maxpv>0, v=pv/maxpv; else v=zeros(nz,1); end

% remove cybernetic variables
if ~isempty(strmatch('u',cyb.vars2fix,'exact'))
    u = ones(nz,1);
end
if ~isempty(strmatch('v',cyb.vars2fix,'exact'))
    v = ones(nz,1);
end
if ~isempty(strmatch('erel',cyb.vars2fix,'exact'))
    erel = ones(nz,1);
end

rz=v.*erel.*rkin;

c=x(iBiom); % biomass

mu=model.sxz(iBiom,:)*rz;
mumax = model.sxz(iBiom,:)'.*kmax;
% mumax_ = model.sxz(iBiom,:)*kmax;
% mumax = mumax_*ones(nz,1);

dydt(1:nx)=model.sxz*rz*c;
dydt(nx+1:nx+nz)=(mumax+be)./(al+ke).*(al+u.*rekin)-diag(be)*erel-mu*erel;
dydt(ny) = flowin - flowout;

dxdt_add = flowin*(xfeed-x) + kla.*(xsat-x) + dxdt_more;
dxdt_add = dxdt_add/vol;

% additional terms
dydt(1:nx)=dydt(1:nx)+dxdt_add;

function [rz,u,v,erel]=getEmvars(t,y,model)
kin = model.kin;
cyb = model.cyb;

nx = length(model.x); ny = length(y); nz = size(model.z,2);
kmax = kin.kmax.num; al=cyb.al;be=cyb.be; ke=cyb.ke;
iBiom = strmatch(model.biom,model.x,'exact');

dydt=zeros(ny,1);
for i=1:ny
    if y(i)<0,y(i)=eps;end
end
x=y(1:nx);
erel=y(nx+1:nx+nz);

% rkin from a user-defined function (basic)
[rkin,flowin,flowout,xfeed,kla,xsat,dxdt_more] = udf_basic(t,kin,x,model);
% rekin
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
        roi=cyb.ncarbon.*roi;
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
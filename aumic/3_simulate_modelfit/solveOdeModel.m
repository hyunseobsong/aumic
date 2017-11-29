function [result,model] = solveOdeModel(model,y0,tspan,kinetic,cybernetic)

kin.kmax.num = cell2mat(kinetic.kmax(:,2));
kin.kmax.name = kinetic.kmax(:,1);
kin.K.num = cell2mat(kinetic.K(:,2));
kin.K.name = kinetic.K(:,1);
cyb = cybernetic;

model.kin = kin;
model.cyb = cyb;

options4ode = [];
[T,Y] = ode15s(@odeModel,tspan,y0,options4ode,kin,cyb,model);

nz = size(model.z,2); nT = length(T); nr = size(model.r,1);
RZ = zeros(nT,nz); U = RZ; V = RZ; EREL = RZ; FLUX = zeros(nT,nr);
for i=1:nT
    t = T(i); y = Y(i,:)';
    [rz,u,v,erel] = getEmvars(t,y,model);
    RZ(i,:) = rz';
    U(i,:) = u';
    V(i,:) = v';
    EREL(i,:) = erel';
    FLUX(i,:) = (model.z*rz)';
end

result.T = T;
result.Y = Y;
result.RZ = RZ;
result.U = U;
result.V = V;
result.EREL = EREL;
result.FLUX = FLUX;

function cybernetic = checkSimulationSetting(model,initial,kinetic,cybernetic)

% check the dimension of kmax
nz = size(model.z,2);
[nrow,ncol]  = size(kinetic.kmax);
if nrow ~= nz || ncol ~=2
    error('The dimension of kmax should be %d-by-2!',nz)
end

% check metaObj and ncarbon
switch cybernetic.metabObj
    case 'carbon uptake'
        cybernetic.ncarbon = cybernetic.ncarbon(:);
        nncarbon  = size(cybernetic.ncarbon,1);
        if nncarbon ~= nz
            error('The length of ''ncarbon'' should be %d!',nz)
        end
    case 'growth'
        if ~isempty(cybernetic.ncarbon)
            warning('''cybernetic.ncarbon'' will be neglected.')
        end
    otherwise
        error ('cybernetic.metaObj ''%s'' is not supported!',cybernetic.metabObj)
end

% check al,be,ke
cybernetic.al = cybernetic.al(:);
cybernetic.be = cybernetic.be(:);
cybernetic.ke = cybernetic.ke(:);

nal = size(cybernetic.al,1);
nbe = size(cybernetic.be,1);
nke = size(cybernetic.ke,1);

if nal==1 && nal<nz
    cybernetic.al = cybernetic.al*ones(nz,1);
elseif nal>1 && nal<nz
	error('The length of ''al'' should be %d or scalar!',nz)
end
if nbe==1 && nbe<nz
    cybernetic.be = cybernetic.be*ones(nz,1);
elseif nbe>1 && nbe<nz
	error('The length of ''be'' should be %d or scalar!',nz)
end
if nke==1 && nke<nz
    cybernetic.ke = cybernetic.ke*ones(nz,1);
elseif nke>1 && nke<nz
	error('The length of ''ke'' should be %d or scalar!',nz)
end

% check initial.met name 
init.met.num = cell2mat(initial.met(:,2));
init.met.name = initial.met(:,1);

for imet = 1:length(init.met.num)
    imet_check = find(strcmp(init.met.name{imet},model.x));
    if isempty(imet_check)
        error('initial.met name ''%s'' is not found!',init.met.name{imet})
    end
end

% check initial.relEnz
initial.relEnz = initial.relEnz(:);
nrelEnz = length(initial.relEnz);
if nrelEnz ~= nz
    error('The length of initial.relEnz should be %d!',nz)
end

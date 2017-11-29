function model = addMets(model,met)

load myModeli

met = cellstr(met);
nmet = length(met);
iimet = [];
for imet = 1:nmet
    imet_ = find(strcmp(met{imet},modeli.x));
    if isempty(imet_)
        error('The metabolite ''%s'' is not available for addition!',met{imet})
    end
    imet__ = find(strcmp(met{imet},model.x));
    if ~isempty(imet__)
        error('The metabolite ''%s'' already exists!',met{imet})
    end
    model.x(end+1) = modeli.x(imet_);
    model.sx(end+1,:) = modeli.sx(imet_,:);
    model.sxz(end+1,:) = model.sx(end,:)*model.z;
end

disp(' ')
disp('*************************************************************************')
disp('Be sure the changed list of metabolites should be reflected in udf_basic!')
disp('*************************************************************************')

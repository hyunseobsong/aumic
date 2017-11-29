function model = removeMets(model,met)

met = cellstr(met);
nmet = length(met);
iimet = [];

for imet = 1:nmet
    imet_ = find(strcmp(met(imet),model.x));
    if isempty(imet_)
        error('The metabolite ''%s'' is not available for removal!',met{imet})
    end
    iimet = [iimet;imet_];
end
model.x(iimet) = [];
model.sx(iimet,:) = [];
model.sxz(iimet,:) = [];

disp(' ')
disp('*************************************************************************')
disp('Be sure the changed list of metabolites should be reflected in udf_basic!')
disp('*************************************************************************')

function model = classifyEms(netmodel,famRule,met)

emRule = famRule(:,1);
refFlux = famRule(:,2);

nfam = size(emRule,1);
nZmas = size(netmodel.z,2);

rr = zeros(nfam,nZmas);

zlogical = logical(netmodel.z);

for iFam = 1:nfam
    rule = emRule{iFam};
	rule = strrep(rule,'AND','and');
	rule = strrep(rule,'OR','or');
	rule = strrep(rule,'NOT','not');

    w = regexp(rule,'\<\w*\>','match'); 
    w = setdiff(w,{'and','or','not'});
    rule = strrep(rule,'and','&');
    rule = strrep(rule,'or','|');
    rule = strrep(rule,'not','~');

    for kk = 1:length(w)
        iircheck = find(strcmp(w{kk},netmodel.r));
        if isempty(iircheck)
            error('Reaction name ''%s'' in famRule is not found!',w{kk})
        else
            iiR(kk,1) = iircheck;
        end
        rule = strrep(rule,w{kk},['z_(',num2str(iiR(kk)),')']);
    end
   
    for iZmas = 1:nZmas
        z_ = zlogical(:,iZmas);
        rr(iFam,iZmas) = eval(rule);
    end   

end 

zfam = [];
for iFam = 1:nfam
    iiZfam = find(rr(iFam,:));
    zfam_ =  netmodel.z(:,iiZfam);
    nzInFam(iFam,1) = size(zfam_,2);
    % normalize ems
    irefluxcheck = find(strcmp(refFlux{iFam},netmodel.r));
    if isempty(irefluxcheck)
        error('Reference reaction ''%s'' in famRule is not found!',refFlux{iFam})
    else
        iRefFlux = irefluxcheck;
    end
    for j = 1:nzInFam(iFam,1)
        zfam__(:,j) = zfam_(:,j)/zfam_(iRefFlux,j);
    end    
    % concatenate families of ems
    zfam = [zfam zfam__];
end

% remove extracellur metabolites with no dynamics
sxz_ = netmodel.sx*zfam;
sxz_sum = sum(abs(sxz_),2);
iiX = find(sxz_sum);
sxz = sxz_(iiX,:);
sx = netmodel.sx(iiX,:);


metbiomcheck = find(strcmp(met.biom,netmodel.x));
if isempty(metbiomcheck)
    error('met.biom name ''%s'' is not found!',met.biom)
else
    metBiom = met.biom;
end

% save model
model.x = netmodel.x(iiX); % extracellular metabolites in ODEs 
model.r = netmodel.r;
model.sx = sx;
model.z = zfam;
model.sxz = sxz;
model.nfam = nfam;
model.nzInFam = nzInFam;
model.biom = metBiom;


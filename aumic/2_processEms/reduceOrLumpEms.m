function [zInFam,model] = reduceOrLumpEms(reduceOrLump,iFam,model,data,options)

% check inputs
if nargin < 4
    options = [];
end

nzInFam = [0;model.nzInFam];
z = model.z(:,sum(nzInFam(1:iFam))+1:sum(nzInFam(1:iFam))+nzInFam(iFam+1));

nx = size(data.fam{iFam},1);
if nx>=1
    dataName = data.fam{iFam}(:,1);
    dataValue = cell2mat(data.fam{iFam}(:,2));
end

for ix = 1:nx
    idatacheck = find(strcmp(dataName{ix},model.r));
    if isempty(idatacheck)
        error('Reaction name ''%s'' in emProces.data is not found!',dataName{ix})
    else
        iData(ix,1) = idatacheck;
    end
end

% remove reference row from yexp and iData
if nx>=1
    iData(1) = []; dataName(1) = []; dataValue(1) = []; 
else
    iData = []; dataName = []; dataValue = []; 
end

switch reduceOrLump
    case 'reduce'
        zInFam = reduceEms(z,iData,dataValue,options);
    case 'lump'
        [zInFam,eta] = lumpEms(z,iData,dataValue,model);
end

if iFam == model.nfam
    disp('EM processing is done!')
else
    disp('Press any key to continue')
    pause
end
    



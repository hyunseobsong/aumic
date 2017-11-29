function model = runEmSoftwaretool(networkFile,emComputationTool)

model = [];

switch emComputationTool
    case 'metatool'
        % metatool 5.1
        % (http://pinguin.biologie.uni-jena.de/bioinformatik/netmodels/...
        %  metatool/metatool5.0/metatool5.0.html)
        ex = metatool(networkFile);
        ex.ems = ex.sub' * ex.rd_ems; % elementary modes
        
    case 'efmtool'
        ex = parse(networkFile);
        mnet.stoich = ex.st;
        mnet.reversibilities = 1-ex.irrev_react;
        % move to efmtool folder (because efmtool should run there)
        rootdir=pwd;
        s = what('efmtool');
        efmtoolpath = s.path;
        cd(efmtoolpath);
        mnet=CalculateFluxModes(mnet);
        ex.ems=mnet.efms;
        disp('--------------------------------')
        nz = size(ex.ems,2);
        fprintf('The number of EMs: %d',nz)
        disp(' ')
        % back to the original folder 
        cd(rootdir);
        
    otherwise
        error('emComputationTool ''%s'' is not supported.',emComputationTool)
end

model.m = ex.int_met; % names of intracellular metabolites 
model.x = ex.ext_met; % names of extracellular metabolites 
model.r = ex.react_name; % reaction names
model.sm = ex.st; % stoichiometric matrix (rows correspond to internal metabolites, columns to reactions)
model.sx = ex.ext; % same structure as st, but rows correspond to external metabolites
model.z = ex.ems; % elementary modes
model.sxz = model.sx*model.z;

function aumic()

% check the file existence
tfiles = {
    'step1_computeEms_template.m'
    'step2_processEms_template.m'
    'step3a_simulate_template.m'
    'step3b_datafit_template.m'
    'step4_plotResults_template.m'
    'udf_basic_template.m'
    };

files = {
    'step1_computeEms.m'
    'step2_processEms.m'
    'step3a_simulate.m'
    'step3b_datafit.m'
    'step4_plotResults.m'
    'udf_basic.m'
    };

refstr = tfiles{1};
filelocation = which(refstr);
irefstr = strfind(filelocation,refstr);
folderlocation = filelocation(1:irefstr-1);

for ifile = 1:length(tfiles)

    if exist(files{ifile},'file')==0
        copyfile([folderlocation,filesep,tfiles{ifile}],files{ifile});
    end
    
end

edit(files{1})
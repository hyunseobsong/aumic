function sys= sbmlModel2metatool(SBMLModel, ext_comparts)
% converts a SBMLModel into a metatool sys struct
% it is assumed that the SBMLModel is free of errors (you can use the
% validate functionality of the SBML toolbox to check this)
% ext_comparts is an optional cell array of compartment ids (strings) that
% can be used to specify if species from the listed compartments should be
% considered external
% for metabolite and reaction names the name attribute is preferably used
% but the id is used when no name is given 
if (SBMLModel.SBML_level == 1)
  id_field_name= 'name';
else
  id_field_name= 'id';
end
specID= {SBMLModel.species(:).(id_field_name)};
specName= {SBMLModel.species(:).name};
no_name= cellfun('isempty', specName);
specName(no_name)= specID(no_name);
specExternal= [SBMLModel.species(:).boundaryCondition] ~= 0;
%A# the 'constant' attribute of species is ignored here because species
%A# that have constant='true' and boundaryCondition='false' are not allowed
%A# as reactants or products in SBML
if any(specExternal)
  disp('The following metabolites are considered external because of boundary conditions:')
  disp(specName(specExternal)');
else
  disp('No explicit external metabolites are defined with boundary conditions.');
end

compartments= {SBMLModel.compartment(:).(id_field_name)};
disp('The following compartments are present in the model:');
disp(compartments');
if nargin < 2
  ext_comparts= false(size(compartments));
  for i= 1:length(compartments)
    key= input(['Consider the metabolites in compartment ', compartments{i}, ' as external? y/[n] '], 's');
    if strcmp(key, 'y')
      ext_comparts(i)= true;
    end
  end
  ext_comparts= compartments(ext_comparts);
else
  missing= setdiff(ext_comparts, compartments);
  if ~isempty(missing)
    disp('Error: the following compartments were specified to contain external metabolites but do not exist:');
    disp(missing);
    sys.err= 1;
    return;
  end
end

%A# add additional external species
if ~isempty(ext_comparts) %A# ismember({...}, {}) leads to an error in octave
  specExternal= specExternal | ismember({SBMLModel.species(:).compartment}, ext_comparts);
end

reacID= {SBMLModel.reaction(:).(id_field_name)};
sys.reac_name= {SBMLModel.reaction(:).name};
no_name= cellfun('isempty', sys.reac_name);
sys.reac_name(no_name)= reacID(no_name);
sys.irrev_react= double([SBMLModel.reaction(:).reversible] == 0); %A# don't use logical

%A# prepare one-stop look-up of all species references
num_cons= cellfun('length', {SBMLModel.reaction(:).reactant});
num_prod= cellfun('length', {SBMLModel.reaction(:).product});
spec= cell(1, sum(num_cons) + sum(num_prod));
begin= 1;
for i= 1:length(SBMLModel.reaction)
  [spec(begin:begin+num_cons(i)-1)]= {SBMLModel.reaction(i).reactant(:).species};
  begin= begin + num_cons(i);
  [spec(begin:begin+num_prod(i)-1)]= {SBMLModel.reaction(i).product(:).species};
  begin= begin + num_prod(i);
end
[spec, spec_ind]= ismember(spec, specID); %A# perform look-up
clear spec;

%A# set up stoichiometric matrix
st= zeros(length(specName), length(sys.reac_name));
has_denominators= ~(SBMLModel.SBML_level == 2 && SBMLModel.SBML_version > 1); %A# future-proof logic???
begin= 1;
for i= 1:length(SBMLModel.reaction)
  reac_range= begin:begin+num_cons(i)-1;
  begin= begin + num_cons(i);
  prod_range= begin:begin+num_prod(i)-1;
  begin= begin + num_prod(i);
  st(spec_ind(reac_range), i)= -double([SBMLModel.reaction(i).reactant(:).stoichiometry]);
  st(spec_ind(prod_range), i)= double([SBMLModel.reaction(i).product(:).stoichiometry]);
  if has_denominators
    if ~isempty(reac_range) %A# prevents empty-matrix-dimension-mismatch Matlobotomy
      st(spec_ind(reac_range), i)= st(spec_ind(reac_range), i) ./ double([SBMLModel.reaction(i).reactant(:).denominator]');
    end
    if ~isempty(prod_range) %A# prevents empty-matrix-dimension-mismatch Matlobotomy
      st(spec_ind(prod_range), i)= st(spec_ind(prod_range), i) ./ double([SBMLModel.reaction(i).product(:).denominator]');
    end
  end
end
sys.st= st(~specExternal, :);
sys.ext= st(specExternal, :);
sys.int_met= specName(~specExternal);
sys.ext_met= specName(specExternal);
sys.err= 0;

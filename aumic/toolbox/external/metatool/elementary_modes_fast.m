%% Copyright (C) 2005 Axel von Kamp
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% implementation of the Wagner alogrithm
% the parameter ems is the kernel, the return value ems is the matrix of elementary modes
% irrev_react(i) is 0 when row i of ems corresponds to a reversible reaction
% req_react is the number of reactions at the end of the tableau that are required
% to participate in all final modes

function [ems, irrev]= elementary_modes_fast(ems, irrev_react, ersatz_rd, subsys_rows, rd, req_react)
  if nargin < 6
    req_react= 0;
  end
  irrev= zeros(1, size(ems, 2)); %# all modes start as reversible
  poss_comb= 0; %# possible combinations
  all_reacts_irrev= all(irrev_react); %# true if all reactions are irreversible
  ems_size= size(ems, 2);
  num_ems= ems_size;
  kn= ems;
  kn_cols= size(ems, 2);
  for c= 1:size(ems, 1)
    if c > kn_cols
      fprintf('Row %d; so far %ld possible combinations and %d modes\n', c, poss_comb, size(ems, 2));
    end

    if c == size(ems, 1) + 1 - req_react
      req_react= req_react - 1;
      fprintf('Keeping only modes that incorporate reaction %d\n', c);
      if irrev_react(c)
	keep= find(ems(c, 1:num_ems) < 0 & ~irrev); % find still reversible modes with an entry < 0
	ems(:, keep)= -ems(:, keep); % and turn them around
	keep= find(ems(c, 1:num_ems) > 0);
	irrev= ones(1, length(keep));
      else
	keep= find(ems(c, 1:num_ems) ~= 0);
	irrev= irrev(1, keep);
      end
      ems= ems(:, keep);
      num_ems= length(keep);
      ems_size= num_ems;
      continue;
    end

    cprev= c - 1;

    %current_subsys= kernel_fp(kn(1:c,:)')';
    %if size(current_subsys, 1) < c
    %  req_zeros= c - size(current_subsys, 1);
    %else
    %  req_zeros= 0;
    %end
    current_subsys= ersatz_rd(1:subsys_rows(c), :);
    if subsys_rows(c) < c
      req_zeros= c - subsys_rows(c);
    else
      req_zeros= 0;
    end

    req_zeros= req_zeros - 2;
    pos_ind= (ems(c, :) > 0);
    neg_ind= (ems(c, :) < 0);
    zero_ind= ~(pos_ind | neg_ind);
    rev_ind= find(~(zero_ind | irrev));
    irrev_pos_ind= find(pos_ind & irrev);
    irrev_neg_ind= find(neg_ind & irrev);
    zero_ind= find(zero_ind);

    % restructure tableau
    irrev_neg_first= (length(irrev_neg_ind) > 0) * 1.0; % make this scalar
    irrev_neg_last= length(irrev_neg_ind);
    irrev_pos_first= irrev_neg_last + (length(irrev_pos_ind) > 0);
    irrev_pos_last= irrev_neg_last + length(irrev_pos_ind);
    rev_first= irrev_pos_last + (length(rev_ind) > 0);
    rev_last= irrev_pos_last + length(rev_ind);
    ems= [ems(:, irrev_neg_ind), ems(:, irrev_pos_ind), ems(:, rev_ind), ems(:, zero_ind)];
    irrev= [ones(1, length(irrev_pos_ind) + length(irrev_neg_ind)), zeros(1, length(rev_ind)), irrev(zero_ind)];

    %# collect combinations with sufficient common zeros
    fprintf('Counting common zeros...');
    comb_size= 32;
    combinations= zeros(2, comb_size); %# collects ems indices
    comb_ind= 0;
    if length(rev_ind) > 1 %# REV/REV combinations
      poss_comb= poss_comb + length(rev_ind) * (length(rev_ind) - 1);
      for i= rev_first:(rev_last - 1)
	for j= (i+1):rev_last
	  if sum(~(ems(1:cprev, i) | ems(1:cprev, j))) >= req_zeros
	    comb_ind= comb_ind + 1;
	    if comb_ind > comb_size %# reallocate if neccesary
	      comb_size= comb_size * 2;
	      combinations(:, comb_ind:comb_size)= 0;
	    end 
	    combinations(:, comb_ind)= [i, j]';
	  end
	end
      end
    end %if
    if length(rev_ind) > 0
      if length(irrev_pos_ind) > 0 %# IRREV_POS/REV combinations
	poss_comb= poss_comb + length(rev_ind) * length(irrev_pos_ind);
	for i= rev_first:rev_last
	  for j= irrev_pos_first:irrev_pos_last
	    if sum(~(ems(1:cprev, i) | ems(1:cprev, j))) >= req_zeros
	      comb_ind= comb_ind + 1;
	      if comb_ind > comb_size %# reallocate if neccesary
		comb_size= comb_size * 2;
		combinations(:, comb_ind:comb_size)= 0;
	      end 
	      combinations(:, comb_ind)= [i, j]';
	    end
	  end
	end
      end %if
      if length(irrev_neg_ind) > 0 %# IRREV_NEG/REV combinations
	poss_comb= poss_comb + length(rev_ind) * length(irrev_neg_ind);
	for i= rev_first:rev_last
	  for j= irrev_neg_first:irrev_neg_last
	    if sum(~(ems(1:cprev, i) | ems(1:cprev, j))) >= req_zeros
	      comb_ind= comb_ind + 1;
	      if comb_ind > comb_size %# reallocate if neccesary
		comb_size= comb_size * 2;
		combinations(:, comb_ind:comb_size)= 0;
	      end 
	      combinations(:, comb_ind)= [i, j]';
	    end
	  end
	end
      end %if
    end %if
    %# IRREV_POS/IRREV_NEG combinations
    if length(irrev_pos_ind) > 0 & length(irrev_neg_ind) > 0 
      poss_comb= poss_comb + length(irrev_neg_ind) * length(irrev_pos_ind);
      for i= irrev_pos_first:irrev_pos_last
	for j= irrev_neg_first:irrev_neg_last
	  if sum(~(ems(1:cprev, i) | ems(1:cprev, j))) >= req_zeros
	    comb_ind= comb_ind + 1;
	    if comb_ind > comb_size %# reallocate if neccesary
	      comb_size= comb_size * 2;
	      combinations(:, comb_ind:comb_size)= 0;
	    end 
	    combinations(:, comb_ind)= [i, j]';
	  end
	end
      end
    end %if
    
    combinations= combinations(:, 1:comb_ind); %# resize
    fprintf(' reducing combination candidates...');
    comb_pattern= ems(1:cprev, combinations(1, :)) | ems(1:cprev, combinations(2, :));
    [comb_pattern, sv, ind]= sort_modes(comb_pattern);
    combinations= combinations(:, ind);
    ind= unique_indices(sv);
    comb_pattern= comb_pattern(:, ind);
    combinations= combinations(:, ind);

    fprintf(' doing rank tests...');
    for i= 1:size(combinations, 2)
      current_pat= comb_pattern(:, i);
      submat= current_subsys(:, find(current_pat));
%      if prod(size(submat)) == 0 %# prevents an error in octave's rank(); error seems to be fixed
%	rg= 0;
%      else
%	rg= rank(submat);
%      end
      if rank(submat) == sum(current_pat) - 1
	r1= combinations(1, i);
	r2= combinations(2, i);
	f1= -ems(c, r1);
	f2= ems(c, r2);
	if (irrev(1, r1) & f1 < 0) | (irrev(1, r2) & f2 < 0)
	  f1= f1 * -1;
	  f2= f2 * -1;
	end
	num_ems= num_ems + 1;
	if num_ems > ems_size %# reallocate if neccesary
	  ems_size= ceil(ems_size * 1.5);
	  ems(:, num_ems:ems_size)= 0;
	  irrev(1, num_ems:ems_size)= 0;
	end
	ems(:, num_ems)= ems(:, r1)/f1 + ems(:, r2)/f2;
        ind= find(abs(ems(:, num_ems)) < 1E-10);
	ems(ind, num_ems)= 0;
	if all_reacts_irrev | all((ems(1:cprev, num_ems) ~= 0) == current_pat)
          irrev(1, num_ems)= irrev(1, r1) | irrev(1, r2);
	  ems(:, num_ems)= ems(:, num_ems) / min(1, max(abs(ems(:, num_ems)))); %# normalize new mode
	else
	  num_ems= num_ems - 1;
	end %if
      end %if
    end %for

    fprintf(' updating tableau\n');
    ems= ems(:, 1:num_ems); %# synchronize sizes
    if isempty(ems) % prevents indexing errors
      irrev= [];
    else
      irrev= irrev(1, 1:num_ems);
    end

    if irrev_react(c)
      if length(rev_ind) > 0 % prevents indexing errors
	rev_ind= rev_first:rev_last; % takes restructuring of ems into account
	irrev(1, rev_ind)= 1; %# change to irreversible
	ind= find(ems(c, rev_ind) < 0); %# extract previously reversible modes...
	rev_ind= rev_ind(ind); %# ...with negative coefficients...
	ems(:, rev_ind)= -ems(:, rev_ind); %# ...and change sign
      end
      if length(irrev_neg_ind) > 0 % prevents indexing errors
	irrev_neg_ind= irrev_neg_first:irrev_neg_last; % takes restructuring of ems into account
	%fprintf('Deleting %d modes\n', length(irrev_neg_ind));
	ems(:, irrev_neg_ind)= []; %# delete modes that violate irreversibility
	irrev(irrev_neg_ind)= [];
	num_ems= num_ems - length(irrev_neg_ind);
	ems_size= ems_size - length(irrev_neg_ind);
      end
    end %if
    if num_ems == 0
      break; % for c= 1:size(ems, 1)
    end
  end % for c= 1:size(ems, 1)
  ems= ems(:, 1:num_ems); %# synchronize sizes
  if isempty(ems) % prevents indexing errors
    irrev= [];
  else
    irrev= irrev(1, 1:num_ems);
  end
  fprintf('Found %d elementary modes; examined %ld possible combinations\n', size(ems, 2), poss_comb);
%endfunction

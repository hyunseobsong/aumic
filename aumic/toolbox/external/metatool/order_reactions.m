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

% do depth first search through the network
% st must contain the entries for the external metabolites that occur in start_metabs
% and end_metabs; make sure the indices refer to the correct rows!
% what happens when start_metabs and end_metabs overlap?

function order= order_reactions(st, irr, start_metabs, end_metabs, met_name, react_name)
  branches= reshape(start_metabs, length(start_metabs), 1); % column vector
  branches= flipud(branches);
  end_met= zeros(1, size(st, 1));
  end_met(end_metabs)= 1;
  %is_ext_met= end_met'; % the metabolites that are regarded as extern in this script
  %is_ext_met(start_metabs)= 1;
  unused_metab= ones(size(st, 1), 1); % column vector
  unused_metab(start_metabs)= 0;
  unused_metab(end_metabs)= 0;
  unused_react= ones(1, size(st, 2)); % row vector
  order= zeros(1, size(st, 2));
  order_ind= 1;
  metab_count= sum(st' ~= 0);
%#  [tmp, ind]= sort(metab_count(branches));
%#  branches= branches(ind(end:-1:1)); % put least often used metabolite last

  while ~isempty(branches)
    metab= branches(length(branches));
    %disp(met_name(metab));
    if end_met(metab)
      branches(length(branches))= [];
      continue;
    end
    cons_react= (st(metab, :) ~= 0) & ((st(metab, :) < 0) | ~irr) & unused_react;
    %disp('Branches:');
    %disp(react_name(find(cons_react)));
    if any(cons_react)
      % follow the reaction that produces the least number of unused metabolites
      cons_react= find(cons_react);
      prod_metabs= zeros(1, length(cons_react));
      for ri= 1:length(cons_react)
	react= cons_react(ri);
	prod_metabs(ri)= sum(sign(st(find(unused_metab), react)) == -sign(st(metab, react)));
      end
      %## folgendes kann man auch durch min ersetzen
      %#[tmp, ind]= sort(prod_metabs);
      %#react= cons_react(ind(1));
      [tmp, ind]= min(prod_metabs);
      react= cons_react(ind);
      %#react= cons_react(1);
      %disp(react_name(react));
      order(order_ind)= react;
      order_ind= order_ind + 1;
      unused_react(react)= 0;
      % append all unused metabolites that are produced by react to branches list
      prod_metab= (sign(st(:, react)) == -sign(st(metab, react))) & unused_metab;
      unused_metab(prod_metab)= 0;
      %branches= [find(prod_metab & is_ext_met); branches; find(prod_metab & ~is_ext_met)];
      prod_metab= find(prod_metab);
      [tmp, ind]= sort(metab_count(prod_metab));
      branches= [branches; prod_metab(fliplr(ind))]; % put least often used metabolite last
    else
      branches(length(branches))= [];
    end
  end
  order=order(find(order)); % not all reactions may have been visited

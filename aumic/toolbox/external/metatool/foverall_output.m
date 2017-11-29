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

% when incl_int is ~= 0 then the internal metabolites are included as well

function foverall_output(fid, fluxmat, sys, incl_int)
  if fid == -1
    return;
  end
  if nargin < 4
    incl_int= 0;
  end
  fprintf(fid, ' overall reaction\n\n');
  if incl_int == 0
    if isempty(sys.ext_met)
      fprintf(fid, ' - not found -\n\n');
      return;
    end
    names= sys.ext_met;
    overall= round(sys.ext * fluxmat * 16384) / 16384; % clumsy rounding
  else
    names= cell(size(sys.int_met, 1) + size(sys.ext_met, 1), 1);
    ind= 1; % clumsy merging follows...
    for i= 1:size(sys.int_met, 1)
      names{ind}= sys.int_met{i};
      ind= ind + 1;
    end
    for i= 1:size(sys.ext_met, 1)
      names{ind}= sys.ext_met{i};
      ind= ind + 1;
    end
    overall= round([sys.st; sys.ext] * fluxmat * 10000) / 10000; % clumsy rounding
  end
  for j= 1:size(overall, 2)
    consumed= find(overall(:, j) < 0)'; % transpose to make a row vector
    produced= find(overall(:, j) > 0)'; % transpose to make a row vector
    fprintf(fid, ' %d:\t', j);
    if (length(consumed) + length(produced)) == 0
      fprintf(fid, 'no net transformation of external metabolites\n');
      continue;
    end
    plus= length(consumed) - 1;
    for m= consumed
      if overall(m, j) ~= -1
	fprintf(fid, '%g ', -overall(m, j));
      end
      fprintf(fid, '%s ', names{m});
      if plus > 0
	fprintf(fid, '+ ');
	plus= plus - 1;
      end
    end
    %if irr(j)
    %  fprintf(fid, '-> ');
    %else
      fprintf(fid, '= ');
    %end
    plus= length(produced) - 1;
    for m= produced
      if overall(m, j) ~= 1
	fprintf(fid, '%g ', overall(m, j));
      end
      fprintf(fid, '%s ', names{m});
      if plus > 0
	fprintf(fid, '+ ');
	plus= plus - 1;
      end
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '\n');
  %fflush(fid);

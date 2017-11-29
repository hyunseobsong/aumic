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

function fenzyme_output(fid, mat, irrev_ems, react_name, offset)
  if fid == -1
    return;
  end
  if nargin < 5
    offset= 0;
  end
  if offset == 0
    fprintf(fid, ' enzymes\n');
  end
  mat= round(mat * 16384) / 16384; % clumsy rounding
  for j= 1:size(mat, 2)
    ind= find(mat(:, j))'; % transpose to make a row vector
    fprintf(fid, '\n %d: (%d)\t', j + offset, length(ind));
    for r= ind
      if abs(mat(r, j)) == 1
	if mat(r, j) < 0
	  fprintf(fid, '-');
	end
      else
	fprintf(fid, '%g ', mat(r, j));
      end
      fprintf(fid, '%s ', react_name{r});
    end
    if irrev_ems(j) == 1
      fprintf(fid, 'irreversible');
    elseif irrev_ems(j) == 0
      fprintf(fid, 'reversible');
    end
  end
  if nargin < 5
    fprintf(fid, '\n\n');
  end
  %fflush(fid);

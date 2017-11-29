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

function fcons_rel_output(fid, crel, int_met)
  if fid == -1
    return;
  end
  crel= round(crel * 16384) / 16384; % clumsy rounding
  
  for j= 1:size(crel, 2)
    fprintf(fid, ' %d:\t', j);
    ind= find(crel(:, j))'; % transpose to make a row vector
    plus= length(ind) - 1;
    for m= ind
      if abs(crel(m, j)) == 1
	if crel(m, j) < 0
	  fprintf(fid, '-');
	end
      else
	fprintf(fid, '%g ', crel(m, j));
      end
      fprintf(fid, '%s ', int_met{m}); % cannot find a less clumsy way...
      if plus > 0
	fprintf(fid, '+ ');
	plus= plus - 1;
      end
    end
    fprintf(fid, '= const\n');
  end
  fprintf(fid, '\n');

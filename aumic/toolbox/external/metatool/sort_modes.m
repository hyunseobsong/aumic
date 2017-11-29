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

function [ems, s, ind] = sort_modes(ems)
  if isempty(ems)
    s= zeros(0, size(ems, 2));
    ind= s;
    return;
  end
  [m, n]= size(ems);
  if m > 52
    error('Cannot sort modes with more than 52 rows\n');
  end
  bin= pow2((1:m) - 1)';
  val= zeros(1, n);
  for j=1:n
    val(j)= sum((ems(:, j) ~= 0) .* bin); % the binary number represented by this mode
  end
  [s, ind]= sort(val);
  ems= ems(:, ind);
%endfunction

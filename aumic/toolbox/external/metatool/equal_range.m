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

% finds equal ranges in sorted vector vec
% the equal ranges begin at the indices in rav (last index is length(vec) + 1)

function rav= equal_range(vec, tol)
  if nargin < 2
    tol= 1e-10;
  end
  if isempty(vec)
    rav= [];
    return;
  end
  rav= [1];
  for i= 1:(length(vec) - 1)
    if abs(vec(i) - vec(i + 1)) > tol
      rav(length(rav) + 1)= i + 1;
    end
  end
  rav(length(rav) + 1)= length(vec) + 1;

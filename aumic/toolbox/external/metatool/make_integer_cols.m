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

function I = make_integer_cols(K)
%#  if nargin < 2
%#    tol = 1e-6 * norm(K, 1)
%#  end
%#  [n, d]= rat(K, tol);
  I= zeros(size(K));
  for i= 1:size(K, 2)
    %#I(:, i)= (n(:, i) / gcd(n(:, i))) .* (lcm(d(:, i)) ./ d(:, i));
    tol= 1e-6 * sum(abs(K(:, i)));
    again= 1;
    while tol > 1e-20 & again
      [ni, di]= rat(K(:, i), tol); %# rat erfüllt die Toleranz komponentenweise
      ind= find(ni);
      if var((ni(ind) ./ di(ind)) ./ K(ind, i)) < 1e-20
	again= 0;
      else
	tol= tol / 2;
      end
    end
    I(:, i)= (ni / fast_gcd(abs(ni))) .* (local_lcm(di) ./ di);

%### die Schleife sollte testen, ob nach Erniedrigung der Toleranz eine bessere Zerlegung vorliegt! Wenn sich die Zerlegung nicht ändert oder schlechter wird, die vorherige Zerlegung als Ergebnis nehmen.

%    while tol > eps & again
%      [ni, di]= rat(K(:, i), tol);
%      I(:, i)= (ni / gcd(ni)) .* (lcm(di) ./ di);
%      ind= find(ni);
%      if var(I(ind, i) ./ K(ind, i)) < eps
%	again= 0;
%      else
%	tol= tol / 2;
%      end
%    end

  end

%## Copyright (C) 1996, 1997 John W. Eaton
%##
%## The following is a part of Octave's lcm function.
%##
%## Octave is free software; you can redistribute it and/or modify it
%## under the terms of the GNU General Public License as published by
%## the Free Software Foundation; either version 2, or (at your option)
%## any later version.
%##
%## Octave is distributed in the hope that it will be useful, but
%## WITHOUT ANY WARRANTY; without even the implied warranty of
%## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%## General Public License for more details.
%##
%## You should have received a copy of the GNU General Public License
%## along with Octave; see the file COPYING.  If not, write to the Free
%## Software Foundation, 59 Temple Place - Suite 330, Boston, MA
%## 02111-1307, USA.

function l = local_lcm (a)
  if (any (a) == 0)
    l = 0;
  else
    a = abs (a);
    l = a (1);
    for k = 1:(length (a) - 1)
      l = l * a(k+1) / gcd (l, a(k+1));
    end
  end

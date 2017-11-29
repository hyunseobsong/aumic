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

function [K, subsys_cols, id_part] = kernel(A)

  [m,n] = size(A);

  [R, pivcol]= rref(A);
  tol= eps * max(m, n) * norm(A, inf); %# same tolerance as used in rref
  r = length(pivcol); % should be rank of A
  if r~= 0 % protect from error in rank
    if r ~= rank(A)
      warning('rref in kernel calculation gives wrong rank, trying workaround...');
      [K, pivcol]= rref(null(A)');
      K= K';
      %#id_part= 1:size(K, 2);
      id_part= pivcol; %# 20.4.2007
      K(abs(K) < tol)= 0;%# 23.4.2007
      subsys_cols= [];
      return;
    end
  end
  id_part = 1:n;
  id_part(pivcol) = [];
  R(abs(R) < tol)= 0;
  K = zeros(n,n-r);
  if n > r
    K(id_part,:) = eye(n-r,n-r);
    if r > 0
      K(pivcol,:) = -R(1:r,id_part);
    end
  end

  %# set up subsys_cols so that K(1:i, 1:subsys_cols(i)) is a kernel of A(:, 1:i)
  isnopiv= ones(1, n);
  isnopiv(pivcol)= 0;
  subsys_cols= cumsum(isnopiv);

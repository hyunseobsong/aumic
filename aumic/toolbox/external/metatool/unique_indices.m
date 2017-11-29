%# sv is a vector of sorted values; sv(ui) are unique values

function ui= unique_indices(sv)
  ui= 1:length(sv);
  src= 1;
  dest= 1;
  while src < length(sv)
    src= src + 1;
    if sv(ui(src)) ~= sv(ui(dest))
      dest= dest + 1;
      ui(dest)= ui(src);
    end %if
  end
  ui((dest+1):length(ui))= [];
%endfunction

## Copyright (C) 2023 Carlo de Falco
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

function [reactions, index] = read_reactions (filename)

  fid = fopen (filename, 'r');
  str = fread (fid, inf, 'char');
  fclose (fid);
  reactions = jsondecode (char (str.'), "makeValidName", false);

  for the_reaction = reactions(:)'
    for [value, key] = the_reaction.reactants
      index.(key) = 0;
    endfor
    for [value, key] = the_reaction.products
      index.(key) = 0;
    endfor
  endfor
  idx = 1;
  for [~, key] = index
    index.(key) = idx++;
  endfor
  
endfunction

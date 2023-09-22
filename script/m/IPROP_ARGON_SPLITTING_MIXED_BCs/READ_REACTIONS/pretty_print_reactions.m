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

function pretty_print_reactions (reactions)
  for the_reaction = reactions(:).'

    for [value, key] = the_reaction.reactants
      if (value > 1) printf ("%3.d", value); else printf ("   "); endif
      printf ("%4.4s +", key)
    endfor
    printf ([char(8), "= "])
    
    for [value, key] = the_reaction.products
      if (value > 1) printf ("%3.d", value); else printf ("   "); endif
      printf ("%5.5s +", key)
    endfor
    printf ([char(8) " "])

    printf (" [ ")
    printf ("%5.5g ", the_reaction.rate_coeffs)
    printf ("]\n")
    
  endfor
endfunction

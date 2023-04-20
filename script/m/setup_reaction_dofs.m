function [dofnames] = setup_reaction_dofs (reactions)
  for the_reaction = reactions(:)'
    for [value, key] = the_reaction.reactants
      s.(key) = 0;
    endfor
    for [value, key] = the_reaction.products
      s.(key) = 0;
    endfor
  endfor
  dofnames = fieldnames (s);
endfunction

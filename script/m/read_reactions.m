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

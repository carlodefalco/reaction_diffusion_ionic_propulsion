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

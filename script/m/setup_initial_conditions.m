function ic = setup_initial_conditions (dofnames, ics)
  ic = zeros (numel (dofnames), 1); 
  for ii = 1 : numel (dofnames)
    if (isfield (ics, dofnames{ii}))
      ic(ii) = ics.(dofnames{ii});
    endif
  endfor
endfunction

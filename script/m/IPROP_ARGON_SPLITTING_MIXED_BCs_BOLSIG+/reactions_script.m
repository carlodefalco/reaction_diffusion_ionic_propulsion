%===============================================================================
% READ THE REACTION FILE INPUT
%===============================================================================
addpath (canonicalize_file_name ("./READ_REACTIONS"));
addpath (canonicalize_file_name ("../../../data"));
[r,idx] = read_reactions(file_in_loadpath ("balcon_et_al_argon_ionization.json"));
%pretty_print_reactions (r);
N_species = numfields(idx);
elements = fieldnames(idx);

save -binary -z reactions.gz

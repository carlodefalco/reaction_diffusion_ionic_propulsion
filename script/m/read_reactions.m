function reactions = read_reactions (filename)

  fid = fopen (filename, 'r');
  str = fread (fid, inf, 'char');
  fclose (fid);
  reactions = jsondecode (char (str.'));

endfunction

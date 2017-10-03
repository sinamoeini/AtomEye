function matrix = T2VoightToFull(vector)
  matrix = [ vector(1) vector(3) 
             vector(3) vector(2) ];
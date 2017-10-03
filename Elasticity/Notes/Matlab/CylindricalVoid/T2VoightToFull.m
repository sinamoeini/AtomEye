function matrix = T2VoightToFull(vector)
  matrix = [ vector(1) vector(6) vector(5) 
             vector(6) vector(2) vector(4) 
             vector(5) vector(4) vector(3) ];
  


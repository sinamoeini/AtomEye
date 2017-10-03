% SphericallyIsotropicUnitVectors.m %

function V = SphericallyIsotropicUnitVectors(N)
  theta = acos(2*rand(1,N)-1);
  phi   = rand(1,N)*2*pi;
  V = [ sin(theta).*cos(phi)
        sin(theta).*sin(phi)
        cos(theta) ];

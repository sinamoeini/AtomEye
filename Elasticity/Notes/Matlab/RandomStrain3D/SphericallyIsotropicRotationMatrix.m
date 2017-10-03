% SphericallyIsotropicRotationMatrix.m %

function R = SphericallyIsotropicRotationMatrix()
  V = SphericallyIsotropicUnitVectors(2);
  V(:,2) = V(:,2) - (V(:,1)'*V(:,2))*V(:,1);
  V(:,2) = V(:,2) / norm( V(:,2) );
  R = [ V(:,1) V(:,2) cross(V(:,1),V(:,2)) ];
  
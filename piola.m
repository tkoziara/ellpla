function P = piola (young, poisson, F)

  % Lame's lambda and mi coefficients
  lambda = young*poisson / ((1.0 + poisson)*(1.0 - 2.0*poisson));
  mi = young / (2.0*(1.0 + poisson));

  % J = det (F)
  J = det(F);

  % calculate Green tensor: E = (F'F - 1) / 2
  E = 0.5*(F'*F-eye(3,3));
  
  % obtain the second PK tensor
  % trough the Saint Venant - Kirchhoff law
  S = 2.0 * mi * E;
  trace = E(1,1) + E(2,2) + E(3,3);
  S(1,1) = S(1,1) + lambda * trace;
  S(2,2) = S(2,2) + lambda * trace;
  S(3,3) = S(3,3) + lambda * trace;

  % now conver S - which is the second PK tensor
  % into the first PK tensor: P = F S
  P = F*S;
end
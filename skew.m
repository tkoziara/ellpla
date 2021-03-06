function B = skew(x)

B = zeros(3);

B(1,2) = -1 * x(3);
B(1,3) = x(2);
B(2,1) = x(3);

B(2,3) = -1 * x(1);
B(3,1) = -1 * x(2);
B(3,2) = x(1);
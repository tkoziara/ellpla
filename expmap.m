function R = expmap (O)
  dp = dot(O,O);
  ln = sqrt(dp);
  W = skew(O);
  if ln < 1E-15
    R = eye(3,3);
  else
    R = eye(3,3) + (sin(ln)/ln)*W + ((1-cos(ln))/dp)*W*W;
  end
end

function W = skew(O)
  W = [0, -O(3), O(2);
       O(3), 0, -O(1);
      -O(2), O(1), 0];
end 
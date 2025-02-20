function A= areasquare(r,z)
A1 = trianglearea(r(1:3),z(1:3));
A2 = trianglearea(r(2:4),z(2:4));
A = A1+ A2;
end

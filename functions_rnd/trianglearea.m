function A = trianglearea(r,z) 
A = lengthline(r(1),r(2),z(1),z(2));
B = lengthline(r(2),r(3),z(2),z(3));
C = lengthline(r(3),r(1),z(3),z(1));
S = (A + B + C) /2;
A = sqrt(S*(S-A)*(S-B)*(S-C));
end

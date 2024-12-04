function [ox, oy] = offset_curve(x,y,w,n)
Nx = length(x);
% normalize n in case it was not
for i = 1:Nx
    n(:,i) = n(:,i)./norm(n(:,i),1);
end
% move points in the normal direction with distance w
xy = [x;y];
xy1 = xy - w.*n;
ox = xy1(1,:); oy = xy1(2,:);
end
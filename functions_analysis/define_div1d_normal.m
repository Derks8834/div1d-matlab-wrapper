function [ox, oy, n] = define_div1d_normal(x,y,w,xcb,ycb,ixp,nwall,wxp_list)
Nx = length(x);
wxp_list(1:2) = [30, 50];
%wxp = ceil(Nx/10);

% directions of the line
dx = diff(xcb);
dy = diff(ycb);
ds = (dy.^2 + dx.^2).^0.5;

% normal direction
for i = 1:length(dx)
n(:,i) = [dy(i)  -dx(i)]/ds(i);
end

% find min distance from wall to start interpolating
inw = zeros(2);
dw(1) = 0;
inw(1,2) = 1;
while dw < abs(w(1))
    inw(1,2) = inw(1,2)+1;
    dw = sum(ds(1:inw(1,2)));
end
% index wall 1
inw(1,1:2) = [1 inw(1,2)];
dw(2) = 0;
inw(2,1) = Nx;
while dw < abs(w(end))
    inw(2,1) = inw(2,1)-1;
    dw = sum(ds(inw(2,1):end));
end
inw(2,1:2) = [inw(2,1) Nx];

% find the point with maximum curvature near Xpoint

% wall direction
% adapt normal 
for i = 1:2
    xp = ixp(i);
    if xp < Nx && xp > 1 % only interp if Xpoint is there
        wxp = wxp_list(i);
    % interpolate normals near X-points (takes out sharp turn)
    if i == 1
    in = [max(1,xp-5*wxp) xp-4*wxp  xp+2*wxp xp+3*wxp];
    %tmps = sign(n(1,in(2)))*round(n(1,in(2)),0);
    %tmps2 = sign(n(2,in(2)))*floor(abs(n(2,in(2))));
%     nint(1,1:4) = [n(1,in(1)), tmps, tmps, n(1,in(4))];
%     nint(2,1:4) = [n(2,in(1)), tmps2, tmps2, n(2,in(4))];
     nint(1,1:4) = [n(1,in(1)), -1, -1, n(1,in(4))];
    nint(2,1:4) = [n(2,in(1)), 1, 1, n(2,in(4))];
    else
    in = [xp-3*wxp xp-2*wxp xp+4*wxp min(xp+5*wxp,Nx)];
    %tmps = sign(n(1,in(3)))*round(n(1,in(3)),0);
    %tmps2 = sign(n(2,in(3)))*round(abs(n(2,in(3))),0);
%     nint(1,1:4) = [n(1,in(1)), tmps, tmps, n(1,in(4))];
%     nint(2,1:4) = [n(2,in(1)), tmps2, tmps2, n(2,in(4))];
        nint(1,1:4) = [n(1,in(1)), -1, -1, n(1,in(4))];
    nint(2,1:4) = [n(2,in(1)), 0, 0, n(2,in(4))];
    end
    % interp xdir
    n(1,in(1):in(4)) = interp1(in,nint(1,:),in(1):in(4));
    % interp ydir
    n(2,in(1):in(4)) = interp1(in,nint(2,:),in(1):in(4));
    

    % the above should work as long as the X points do not turn more than
    % 120 deg.

    % change normal to align with target near wall
    if i == 1
      n2w = [nwall(:,1) n(:,inw(i,1))];
    elseif i ==2    
      n2w= [n(:,inw(i,1)) nwall(:,2)];   
    end

    % interp xdir
    n(1,inw(i,1):inw(i,2)) = interp1(inw(i,:),n2w(1,:),inw(i,1):inw(i,2));
    % interp ydir
    n(2,inw(i,1):inw(i,2)) = interp1(inw(i,:),n2w(2,:),inw(i,1):inw(i,2));
    end
end

% normalize n
for i = 1:length(n)
 n(:,i) = n(:,i)./norm(n(:,i),1);
end

% move points in the normal direction    
xy = [x;y];
xy1 = xy - w.*n;
ox = xy1(1,:); oy = xy1(2,:);
end
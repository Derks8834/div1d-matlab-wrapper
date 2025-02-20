function [position] = get_temperature_front(o,grid,Tfront)
% [position] = get_temperature_front(o,in,Tfront) Summary of this function goes here
%  obtain timeseries for temperature front

pause(0.5)
disp('1');
x = grid.x;
Nx = length(x);
[Nt, Nx] = size(o);
position = zeros(1,Nt);
L = max(x);
for i = 1:Nt
position(i) = interpn(x,o(i,:),Tfront,'linear',L)
end
position = 1;
end


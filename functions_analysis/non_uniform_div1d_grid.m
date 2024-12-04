
function [x,xcb] = non_uniform_div1d_grid(L,Nx,dxmin)
      % set-up non equidistant grid for div1d
      for i = 0: Nx
         xnorm(i) =  (i)/ (Nx);
      end
      % calculate the cell boundaries in the grid
      xcb = L * ( (2.0d+0-dxmin)*xnorm - (1.0d+0-dxmin)*xnorm*xnorm);
      % calculate the grid cell centres
      x = (xcb(0:Nx-1) + xcb(1:Nx)) / 2.0d+0;
end

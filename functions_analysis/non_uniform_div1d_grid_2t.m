
function [x,xcb] = non_uniform_div1d_grid_2t(L,Nx,dxmin)
      % set-up non-equidistant grid with two targets: grid is symmetric
      % first define a normalized array running from 0 to 1 at the cell boundaries
      % note that Nx must be even for this to work correctly!
      if logical(rem(Nx,2)) % remainder should be 1 or 0
        error('grid must be even to be symmetric')
      end
      for i = 1: Nx/2+1
         xnorm(i) =  (i-1)/ (Nx/2);
      end
      % calculate the cell boundaries on the right part of the grid grid
      xcb(Nx/2+1:Nx+1) = (L/2.0d+0) + (L/2.0d+0) .* ( (2.0d+0-dxmin)*xnorm(1:Nx/2+1) - (1.0d+0-dxmin).*xnorm(1:Nx/2+1).*xnorm(1:Nx/2+1) );
      % calculate the cell boundaries on the left part of the grid grid
      xcb(1:Nx/2+1)  = (L/2.0d+0) - (L/2.0d+0) .* ( (2.0d+0-dxmin)*xnorm(Nx/2+1:-1:1) - (1.0d+0-dxmin).*xnorm(Nx/2+1:-1:1).*xnorm(Nx/2+1:-1:1) );
      % calculate the grid cell centres
      x = cb2cc(xcb);
end 
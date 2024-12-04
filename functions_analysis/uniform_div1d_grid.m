function [x,xcb] = uniform_div1d_grid(L,Nx)
      % set-up equidistant grid for div1d
      dx = L/Nx;
      delta_x = dx;
      delta_xcb = dx;
      for i = 1:Nx
         x(i) = dx/2.0d+0 +  (i-1)*dx;
         xcb(i) =  (i)*dx;
      end 
      xcb(0) = 0.0d+0;
      return
end 
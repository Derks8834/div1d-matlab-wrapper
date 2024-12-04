function [v_cb] = cc2cb(Nx,v_cc,delta_x, delta_xcb)
%CC2CB interpolates and extrapolates values from cell centers to boundaries
% interpolate [v_cb] = cc2cb(Nx,v_cc,delta_x, delta_xcb)
	v_cb(2:Nx) = ( v_cc(1:Nx-1).*delta_xcb(2:Nx) + v_cc(2:Nx).*delta_xcb(1:Nx-1) ) ./ (2.0d+0*delta_x(1:Nx-1));
       % extrapolate	
	v_cb(1)  = v_cc(1)  -  (v_cc(2)  - v_cc(1)   )/delta_x(1)  * 0.5*delta_xcb(1);
	v_cb(Nx+1) = v_cc(Nx) +  (v_cc(Nx) - v_cc(Nx-1))/delta_x(Nx-1) * 0.5*delta_xcb(Nx); 
end


function plot_div1d(output,input,varargin)
% plot_div1d_output(div1doutput,div1dinput,varargin)
% div1doutput and input and proc obtained by proces_div1d_output.m

% Author: Gijs Derks
% E-mail: g.l.derks@differ.nl
% July 2023

switch input.version{1}
    case 'v4.0.1'
    plotdiv1d_v401(output,input);

    case 'v4.0.2'
    plotdiv1d_v402(output,input);
    
    return
    otherwise
    disp('older versions supported in script below');
    plotdiv1d_below_v4(output,input);
end

end

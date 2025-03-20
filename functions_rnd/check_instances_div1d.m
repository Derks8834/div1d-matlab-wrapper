function [] = check_instances_div1d(varargin)
%CHECK_INSTANCES_DIV1D running from terminal
% if more than "8" are found, wait untill they finish
% give different number through varargin
num_instances = 8;
try num_instances = varargin{1}; catch; end

[~,tmp] = system('pgrep div1d.exe');
Nact = split(tmp);
while length(Nact) > 8 % wait to continue after less then 8 are present
    [~,tmp] = system('pgrep div1d.exe');

    pause(10);
    Nact = split(tmp);
end

end


function [v_cc] = cb2cc(v_cb)
%CB2CC interpolates  values from cell boundaries to centers
% interpolate [v_cc] = cb2cc(v_cb)
    for i = 1:length(v_cb)-1
        v_cc(i) = mean([v_cb(i),v_cb(i+1)]);
    end
end


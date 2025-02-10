function p = minmax(arr)
    p = zeros(1,3);
    p(1) = mean(arr);
    p(2) = min(arr);
    p(3) = max(arr);
end

    
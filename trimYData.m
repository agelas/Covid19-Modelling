function [yend, trimmedY] = trimYData(y)
    yend = y(end, :);
    y(end, :) = [];
    trimmedY = y;
end
% given a y series data in the form of a matrix, separate the matrix into 2
% parts, the y values at the end of matrix (time interval), and the values
% that preceed it. This is done to easily provide initial values for
% segmented model fitting and eventual appended data 
% (if data generated at each segment are directly appended together, it
% would cause duplication of data at the days of the markers

function [yend, trimmedY] = trimYData(y)
    yend = y(end, :);
    y(end, :) = [];
    trimmedY = y;
end
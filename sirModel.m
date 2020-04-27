% function that outputs the y solutions along a time interval given the
% initial y0, the length of the time interval, and the required
% coefficients in the ODE

function y1 = sirModel(Y0, maxTime, b, s)
    tspan = [1:maxTime];
    R0 = b; %How many people each infected person infects
    recv = s; %Recovery rate
    [~, y1] = sirODE(tspan, Y0, R0, recv);
end
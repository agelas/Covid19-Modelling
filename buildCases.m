function [days, cases] = buildCases(growth)
    days = size(growth);
    days = days(1);
    cases = zeros(days,1);
    cases(1) = growth(1,2);
    for i = 2:days
        cases(i) = cases(i-1) + growth(i,2);
    end
end
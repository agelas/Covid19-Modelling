% custom function used to generate the number of cases per day from the
% reported NEW cases per day.
% done by adding new cases number to the number of total cases from the
% previous day

% input growth matrix is in two columns, first column the day number,
% second column is the new cases each day

function [days, cases] = buildCases(growth)
    days = size(growth);
    days = days(1);
    cases = zeros(days,1);
    cases(1) = growth(1,2);
    for i = 2:days
        cases(i) = cases(i-1) + growth(i,2);
    end
end
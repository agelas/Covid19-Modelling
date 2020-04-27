% This script uses the generated segmented model data from launchModel.m
% for each county and plots them altogether for visual comparison
% the y values generated are in percentages, to allow direct comparison
% between counties

% plot casts every county's timeline to the NYC timeline, as NYC has had
% the earliest outbreaks
% as of 4/12, NYC has had covid cases for 42 days
maxsize = 42;

dates = 1:maxsize;
y = dutchessY;
% for each y set, check if the timeline is the same length as the nyc
% timeline
% if not, append a number of 0 entries to the beginning of the county
% timeline in order to make it the same length as the nyc timeline
if (max(size(y)) ~= maxsize)
    prepend = zeros(maxsize - max(size(y)), 1);
    y = [prepend; y];
end
hold on
% then plot on the same plot for visual comparison
plot(dates, y);


% repeat for each county
y = rocklandY;
if (max(size(y)) ~= maxsize)
    prepend = zeros(maxsize - max(size(y)), 1);
    y = [prepend; y];
end
hold on
plot(dates, y);

y = suffolkY;
if (max(size(y)) ~= maxsize)
    prepend = zeros(maxsize - max(size(y)), 1);
    y = [prepend; y];
end
hold on
plot(dates, y);

y = wcY;
if (max(size(y)) ~= maxsize)
    prepend = zeros(maxsize - max(size(y)), 1);
    y = [prepend; y];
end
hold on
plot(dates, y);

y = nycY;
if (max(size(y)) ~= maxsize)
    prepend = zeros(maxsize - max(size(y)), 1);
    y = [prepend; y];
end
hold on
plot(dates, y);

% add legend and axis labels
legend({'Dutchess, Rural', 'Rockland, Suburban', 'Suffolk, Suburban', 'Westchester, Suburban', 'NYC, Urban'})
ylabel('Percentage of county population infected');
xlabel('Days since first NY case');

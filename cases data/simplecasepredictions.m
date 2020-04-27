% this script takes from the original launchModel script alot, but extends
% it by putting together y data for predictions of case growth if the state
% stopped at certain interventions.

%% build cases

% the growth data needs to be imported from the csv files for county
% for this case: NYC

% setting up the number of cases per day
% custom function to build the cases from the day by day growth data
% extract the time span and the number of infected cases each day from the
% custom buildCases function

[days, cases] = buildCases(growth);

%set up initial conditions on day 1
% pop value to be modified for each county
Pop = 8399000;
cases = cases/Pop;
initCoef = 1;
Infec = cases(1)*initCoef;
yinit = [1, Infec, 0];
% at the beginning there are 100% susceptible
% read the first value from the case growth and divide by the total
% population to get a percentage value (Infec)
% 0% recovered


% the processing and model building process uses percentage instead of pop
% numbers, so that the data can be easily extracted to compare between
% counties with different population counts.



%% setup possible R0 and recv

% set up vectors of coefficient values that we want to test (R0, recvs)
numR0 = 100;
numRecv = 20;
minR0 = 0.3;
minRecv = 0.2;
maxR0 = 1;
maxRecv = 0.5;
R0s = linspace(minR0, maxR0, numR0);
recvs = linspace(minRecv, maxRecv, numRecv);
%%

% set the markers in which intervention methods show effect (account for
% incubation period)
% the mark number represents how many days since the first case in the
% country has been identified did the certain state-wide intervention get
% implemented + a 7 day incubation period to find the day in which each
% measure's effects are to take place.

% mark1 - Cuomo urges people who can work from home to do so
% mark2 - Large gatherings are banned (500+ people)
% mark3 - Gatherings are limited to 50, restaurants and other rec
% facilities such as theaters, gyms, bars, etc. are to be closed
% mark4 - Non essential businesses are closed
% mark5 - requires mask wearing in public
mark1 = 14;
mark2 = mark1+4;
mark3 = mark1+8;
mark4 = mark1+12;
mark5 = mark1+24;

%%

% for each identfied time segment, these are the inputs to the fitModel
% function: 
% the R0 and recv coefficient values that we want to test
% the initial value (the final value from the previous segment or the day 1 values)
% the time interval
% the actual data to fit to

% function will output the combination of coefficents (r0 and recv) that
% yields the curve with the best r value, that r value, and the resulting y
% values (all three of SIR).

% store those values separately, and append a trimmed version of the y
% values for each segment into a complete matrix "modelY", which will be
% the complete generated model.

%fitModel and trimYData both custom functions, see respective files
[mark1minR, mark1R0, mark1bestRecv, mark1besty] = fitModel(R0s, recvs, yinit, 1, mark1, cases);
[ymark1, mark1besty] = trimYData(mark1besty);
[mark2minR, mark2R0, mark2bestRecv, mark2besty] = fitModel(R0s, recvs, ymark1, mark1, mark2, cases);
[ymark2, mark2besty] = trimYData(mark2besty);
[mark3minR, mark3R0, mark3bestRecv, mark3besty] = fitModel(R0s, recvs, ymark2, mark2, mark3, cases);
[ymark3, mark3besty] = trimYData(mark3besty);
[mark4minR, mark4R0, mark4bestRecv, mark4besty] = fitModel(R0s, recvs, ymark3, mark3, mark4, cases);
[ymark4, mark4besty] = trimYData(mark4besty);
[mark5minR, mark5R0, mark5bestRecv, mark5besty] = fitModel(R0s, recvs, ymark4, mark4, mark5, cases);
[ymark5, mark5besty] = trimYData(mark5besty);
[mark6minR, mark6R0, mark6bestRecv, mark6besty] = fitModel(R0s, recvs, ymark5, mark5, days, cases);
[endy, mark6besty] = trimYData(mark6besty);

modelY = [mark1besty;mark2besty;mark3besty;mark4besty;mark5besty;mark6besty;endy];

tspan = 1:days;
hold on;

% ^^ everything above is the same as launchModel, and it generates a fitted
% model.

%% no intervention whatsoever

% takes the coefficient values from the beginning of the model fit and runs
% the remainder of the timeline using the same coefficients
% Absolute worst case scenario
[~, ~, ~, stage1y] = fitModel(mark1R0, mark1bestRecv, yinit, 1, days, cases);
y1 = stage1y;
plot(tspan, stage1y(:,2)*Pop);

ylim auto;
xlabel('Days');
ylabel('Infected Cases');



%% no second
% for all other ones, take the original fitted data, cut off at a certain
% mark, and use the last generated coefficient values to simulate the
% remainder of the timeline, repeat for each intervention

[~, ~, ~, stage2y] = fitModel(mark2R0, mark2bestRecv, ymark1, mark1, days, cases);
y2 = [mark1besty;stage2y];
plot(tspan, y2(:,2)*Pop);


%% no third
[~, ~, ~, stage3y] = fitModel(mark3R0, mark3bestRecv, ymark2, mark2, days, cases);
y3 = [mark1besty;mark2besty;stage3y];
plot(tspan, y3(:,2)*Pop);


%% no fourth
[~, ~, ~, stage4y] = fitModel(mark4R0, mark4bestRecv, ymark3, mark3, days, cases);
y4 = [mark1besty;mark2besty;mark3besty;stage4y];
plot(tspan, y4(:,2)*Pop);


%% no fifth
[~, ~, ~, stage5y] = fitModel(mark5R0, mark5bestRecv, ymark4, mark4, days, cases);
y5 = [mark1besty;mark2besty;mark3besty;mark4besty;stage5y];
plot(tspan, y5(:,2)*Pop);


%%
% plot the actual case number (current timeline) for visual comparison
plot(tspan, cases*Pop);

title('COVID Case Prediction With/Without Interventions');
legend({'No Intervention Whatsoever', 'No Intervention Past Large Gathering Ban', 'No Intervention Past Rec Facilities and Restaurants Ban', 'No Intervention Past Non-Essential Business Ban','Current Timeline'});

%% build cases


% setting up the number of cases per day
% i had a fundamnetal misunderstanding of the ode solvers in MATLAB,
% so the data still needs to be in total cases per day
% custom function to build the cases from the day by day growth data
% compiled

[days, cases] = buildCases(growth);

%set up initial conditions on day 1
Pop = 98320;
cases = cases/Pop;
initCoef = 1;
Infec = cases(1)*initCoef;
yinit = [1, Infec, 0];



%% setup possible R0 and recv

% set up vectors of coefficient values that we want to test (R0, recvs)
numR0 = 100;
numRecv = 20;
minR0 = 0.3;
minRecv = 0.2;
maxR0 = 0.8;
maxRecv = 0.5;
R0s = linspace(minR0, maxR0, numR0);
recvs = linspace(minRecv, maxRecv, numRecv);
%%

% set the markers in which intervention methods show effect (account for
% incubation period), random values rn
mark1 = 0;
mark2 = mark1+4;
mark3 = mark1+8;
mark4 = mark1+12;
mark5 = mark1+24;
mark1 = 1;
%%
% this command will run the model fitting for the entirety of the data as
% one segment
%[OVminR, OVbestR0, OVbestRecv, OVbesty] = fitModel(R0s, recvs, yinit, 1, days, cases);

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
%[mark1minR, mark1R0, mark1bestRecv, mark1besty] = fitModel(R0s, recvs, yinit, 1, mark1, cases);
%[ymark1, mark1besty] = trimYData(mark1besty);
mark1minR = 0;
mark1R0 = 0;
mark1bestRecv = 0;
%mark1besty = [0 0 0];
ymark1 = yinit;
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

modelY = [mark2besty;mark3besty;mark4besty;mark5besty;mark6besty;endy];
coefValues = [mark1R0, mark1bestRecv, mark1minR; mark2R0, mark2bestRecv, mark2minR; mark3R0, mark3bestRecv, mark3minR; mark4R0, mark4bestRecv, mark4minR; mark5R0, mark5bestRecv, mark5minR; mark6R0, mark6bestRecv, mark6minR];

%%

%plot to visualize values as a sanity check
tspan = 1:days;
hold on;
plot(tspan, cases*Pop);
%plot(tspan, OVbesty(:,2)*Pop);
plot(tspan, modelY(:,2)*Pop);
ylim auto;
legend('Putnam Data', 'SegmentedModel');
xlabel('Days');
ylabel('Infected Cases');
%legend('Westchester Data', 'ContinuousModel', 'SegmentedModel');
disp(coefValues); 

%%
putnamY = modelY(:,2);
putnamT = tspan;

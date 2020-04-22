%% build cases

[days, cases] = buildCases(growth);
Pop = 967506;
cases = cases/Pop;
initCoef = 20;
Infec = cases(1)*initCoef;
yinit = [1, Infec, 0];

%% setup possible R0 and recv
numR0 = 100;
numRecv = 20;
minR0 = 0.3;
minRecv = 0.2;
maxR0 = 0.8;
maxRecv = 0.5;
R0s = linspace(minR0, maxR0, numR0);
recvs = linspace(minRecv, maxRecv, numRecv);
%%
mark1 = 20;
mark2 = 27;
mark3 = 35;

%%
[OVminR, OVbestR0, OVbestRecv, OVbesty] = fitModel(R0s, recvs, yinit, 1, days, cases);

%%
[mark1minR, mark1R0, mark1bestRecv, mark1besty] = fitModel(R0s, recvs, yinit, 1, mark1, cases);
[ymark1, mark1besty] = trimYData(mark1besty);
[mark2minR, mark2R0, mark2bestRecv, mark2besty] = fitModel(R0s, recvs, ymark1, mark1, mark2, cases);
[ymark2, mark2besty] = trimYData(mark2besty);
[mark3minR, mark3R0, mark3bestRecv, mark3besty] = fitModel(R0s, recvs, ymark2, mark2, mark3, cases);
[ymark3, mark3besty] = trimYData(mark3besty);
[mark4minR, mark4R0, mark4bestRecv, mark4besty] = fitModel(R0s, recvs, ymark3, mark3, days, cases);
[endy, mark4besty] = trimYData(mark4besty);

modelY = [mark1besty;mark2besty;mark3besty;mark4besty;endy];
coefValues = [mark1R0, mark1bestRecv, mark1minR; mark2R0, mark2bestRecv, mark2minR; mark3R0, mark3bestRecv, mark3minR; mark4R0, mark4bestRecv, mark4minR];

%%
tspan = 1:days;
plot(tspan, OVbesty(:,2)*Pop);
ylim auto;
hold on;
plot(tspan, modelY(:,2)*Pop);
plot(tspan, cases*Pop);
legend('OVmodel', 'SegmentedModel', 'data');
disp(coefValues); 

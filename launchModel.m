%1: 100% of the population 'healthy'
%1.03e-6: percent of population that's infected
%0: Recovered population
%50: Number of days, weeks, idek, that the simulation runs


%% build cases

days = size(growth);
days = days(1);
cases = zeros(days,1);
cases(1) = growth(1,2);
for i = 2:days
    cases(i) = cases(i-1) + growth(i,2);
end
Pop = 967506;
cases = cases/Pop;
Infec = cases(1)*250;
Rec = 0;

%% setup possible R0 and recv
numR0 = 1000;
numRecv = 20;
minR0 = 0.3;
minRecv = 0.2;
maxR0 = 0.6;
maxRecv = 0.5;
R0s = linspace(minR0, maxR0, numR0);
recvs = linspace(minRecv, maxRecv, numRecv);
%%
minR = 0;
bestR0 = 0;
bestRecv = 0;
besty = zeros(days, 1);

%%
for i = 1:numR0
    for j = 1:numRecv
        y = sirModel(1,Infec,Rec,days,R0s(i),recvs(j));
        r = sum((cases-y(:,2)).^2);
        if (minR == 0) || (r < minR)
            bestR0 = R0s(i);
            bestRecv = recvs(j);
            minR = r;
            besty = y(:,2);
        end
    end
end

tspan = 1:days;
plot(tspan, besty*Pop);
hold on;
plot(tspan, cases*Pop);
ylim auto;
legend('model', 'data');
disp({minR, bestR0, bestRecv});


%sirModel(967506, 1, 0, 41,  
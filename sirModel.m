function y1 = sirModel(Pop, Infec, Recovd, maxTime, b, s)
tspan = [1:maxTime];
Y0 = [Pop, Infec, Recovd];
R0 = b; %How many people each infected person infects
recv = s; %Recovery rate, this is a completely random number I chose
[~, y1] = sirODE(tspan, Y0, R0, recv);

%y1 = y1.*Pop; %adjusted to population of Westchester
% 
% plot(tspan, y1);
% title('SIR Model');
% legend({'Suceptible', 'Infected', 'Removed'});
% xlabel('Time (days maybe idk honestly)');

end
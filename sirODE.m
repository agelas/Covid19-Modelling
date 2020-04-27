% feeds the input time span, y0, and coefficients into the ode45 solver
% outputs the time span as well as the y values along the time span

function [t, y1] = sirODE(tspan, y0, R0, recv)

    [t, y1] = ode45(@odefun, tspan, y0);
    
    % define the SIR ode functions in the odefun function
    function dydt = odefun(~,y)
        dydt = zeros(3,1);
        
        dydt(1) = -R0 * y(1) * y(2); %dS/dt, susceptible rate
        dydt(2) = R0 * y(1) * y(2) - recv * y(2); %dI/dt, infected rate
        dydt(3) = recv * y(2); %dR/dt, recovered rate
        
        global odeArray 

        odeArray = dydt;
    end


end
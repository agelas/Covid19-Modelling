# probstats
Prob Stats COVID project

launchModel.m: Ideally the only function that is interacted with at the end
of the day when it comes to controlling sirModel.m and sirODE.m.

sirModel.m: Sets up and passes in the initial conditions necessary for
sirODE.m to do its job. Arguments are susceptible % of population, 
infected %, recovered %, time duration for simulation, R0 (how many people
an infected person infects), and recovery rate. Called from launchModel.m.

sirODE.m: Uses ode45 (nonstiff, medium accuracy) to integrate the three 
differential equations responsible that are part of the SIR model. This 
function should never be called by itself and is only called through 
sirModel.m

sirVisualClassDef.m: Creates an object of type sirVisualClassDef. Follows 
the format of obj = sirVisualClassDef(N,KE). N is the number of people in
simulation environment. Make sure it's a square number like 64 or 100. Don't 
worry about KE for now, just keep it to a number between 0.5 and 2 and you 
should be fine.

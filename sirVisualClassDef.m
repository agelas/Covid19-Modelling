
classdef sirVisualClassDef < handle
    properties
        N;          %number of people
        L;          %size of simulation box
        t;          %current number of time steps completed
        dt;         %time step size
        otime;      %number of time steps between outputs
        rc=1;       %cutoff distance for interaction
        pos;        %positions of people in Nx2 array
        vel;        %velocities of people in Nx2 array
        phand;      %list of hadles to all circles representing people
        timeser;    %list of times at which data is recorded
        Eser;       %the energy per person 
        KEser;      %the kinetic energy per person
        PEser;      %potential energy per person
        Infected;   %how many infected people
        infCount;   %how many infected people at the timestep
        chanceInfection; %individual person's chance of getting infected
        logicalInfected; %logical array detailing whether person is infected
        Susceptible;%how many susceptible people
        ACTIVATE_DISTANCING = 0; %if social distancing is taking place in simulation
    end
    
    methods
        %constructor
        function obj=sirVisualClassDef(N,KE,SD) %KE should be between 0.1-0.8 so like 0.5 probably
            obj.N=N;
            obj.L=sqrt(N);
            avg=KE*2;
            obj.t=0;
            obj.dt=0.01;                        %setting the time step
            [posx,posy]=meshgrid([1:sqrt(N)],[1:sqrt(N)]);  %meshgrid for people indices
            posx=reshape(posx,1,N);            %reformatting x positions 
            posy=reshape(posy,1,N);            %reformatting y positions
            obj.pos=[posx',posy'];              %filling in pos array
            obj.pos=obj.pos-0.5;                %shifting indexes to not fall on whole numbers
            obj.vel=randn(obj.N,2)*sqrt(avg);   %distributing random velocities across the people
            obj.infCount = 0;
            obj.chanceInfection = zeros(1, obj.N);
            obj.logicalInfected = zeros(1, obj.N);
            obj.ACTIVATE_DISTANCING=SD;        %if social distancing is taken into account
            
            %Inserting an infected person into the population
            randomUnluckyPerson = floor((obj.N - 1)*rand(1) + 1)
            obj.logicalInfected(randomUnluckyPerson) = 1;
            
        end
        
        function obj=step(obj)
    
            velocityFirst=obj.vel;              %filling in velocityFirst for movement calculations
            
            for i=1:obj.N                       %looping through all other people in simulation
                
                obj.pos(i,1)=normalUpdate(obj,obj.pos(i,1),obj.vel(i,1));  %updating x coordinate
                obj.pos(i,2)=normalUpdate(obj,obj.pos(i,2),obj.vel(i,2));  %updating y coordinate
                
                rij=zeros(1, obj.N);  %initialize radius to atoms within boundary contributing to force
                rijx=zeros(1, obj.N); %initialize x-distance to atoms within boundary
                rijy=zeros(1, obj.N); %initialize y-distance to atoms within boundary
                
                Fix=zeros(1, obj.N);
                Fiy=zeros(1, obj.N);
                
                potentialEnormal=zeros(1, obj.N); %potential energy from within boundaries
                
                checkBoundaryConditions(obj, obj.pos(i,1), obj.pos(i,2),i); %Boundary corrections for people
                
                for j=1:obj.N %Looping through all other people for distance calcs
                       
                    rijx(j)= obj.pos(j,1)-obj.pos(i,1); %x distance is recorded
                    rijy(j)= obj.pos(j,2)-obj.pos(i,2); %y distance is recorded
                    
                    if rijx(j)>obj.L/2  %x contribution through boundary
                        rijx(j)=rijx(j)-obj.L;
                    elseif rijx(j)<-obj.L/2  %x contribution through boundary
                        rijx(j)=rijx(j)+obj.L;    
                    end
                    
                    if rijy(j)>obj.L/2  %y contribution through boundary
                        rijy(j)=rijy(j)-obj.L; 
                    elseif rijy(j)<-obj.L/2  %y contribution through boundary 
                        rijy(j)=rijy(j)+obj.L;
                    end
                    
                    if sqrt((rijx(j))^2+(rijy(j)^2))<=obj.rc   %Checking if distance is within cutoff distance
                        rij(j)=sqrt((rijx(j))^2+(rijy(j)^2));  %If within range, distance is recorded
                        if(obj.logicalInfected(j) == 1) %Person in radius needs to also be infected
                            chanceInf = 0.01;
                            if obj.ACTIVATE_DISTANCING == 1
                                chanceInf = 0.002;
                            end
                            obj.chanceInfection(i) = obj.chanceInfection(i) + chanceInf; 
                        end
                    else
                        rij(j)=0;
                    end
                    
                    if rij(j)==0    %if radius is 0, x and y distances are set to 0 to discard atom
                        Fix(j)=0;
                        Fiy(j)=0;
                    else
                        Fix(j)=24*((2*rij(j))^-14-(rij(j))^-8)*rijx(j); %Potential energy calculation in x
                        Fiy(j)=24*((2*rij(j))^-14-(rij(j))^-8)*rijy(j); %Potential energy calculation in y
                    end    
                end
                
                r = binornd(1,obj.chanceInfection(i)); %Bernoulli distribution is just binomial distribution with N = 1.
                if (r == 1)
                    obj.logicalInfected(i) = 1; %Person is now infected
                end
                
                %finding total potential energy for the i person
                for k=1:length(rij)
                    if rij(k)==0
                        potentialEnormal(k)=0;
                    else
                        potentialEnormal(k)=4*(rij(k).^-12-rij(k).^-6-obj.rc^-12+obj.rc^-6);
                    end
                end
                
                PE = zeros(1, obj.N);
                PE(i)=sum(potentialEnormal)/36; %summing all contributions for total potential energy

                Fix=sum(Fix);   %summing all the x forces
                Fiy=sum(Fiy);   %summing all the y forces

                obj.vel(i,1)=velocityUpdate(obj,obj.vel(i,1),Fix);   %updating x velocity
                obj.vel(i,2)=velocityUpdate(obj,obj.vel(i,2),Fiy);   %updating y velocity

                obj.pos(i,1)=normalUpdate(obj,obj.pos(i,1),obj.vel(i,1));  %updating x coordinate again to make more lively
                obj.pos(i,2)=normalUpdate(obj,obj.pos(i,2),obj.vel(i,2));  %updating y coordinate
                
            end
            velocitySecond=obj.vel; %second velocity to hold updated velocity for kinetic evergy calculations
            velocityMean=(velocityFirst+velocitySecond)./2; %taking the mean of the velocities
            netVelocity=sqrt(velocityMean(:,1).^2+velocityMean(:,2).^2);   %taking the net velocity 
            
            KE=0.5.*(netVelocity.^2);   %calculatingn KE accoring to 0.5mv^2
            KE=KE';                     %transposing KE so it adds properly
            
            Energy(i)=PE(i)+KE(i);      %Calculating total energy
            
            obj.KEser(end+1)=mean(mean(KE)); %filling in KEser array
            obj.PEser(end+1)=mean(PE);       %filling in PEser array
            obj.Eser(end+1)=mean(Energy);    %filling in Eser array
            obj.timeser(end+1)=obj.t;        %filling in timeser array
            obj.Infected(end+1)=obj.infCount;%updating infected count
            obj.Susceptible(end+1)=obj.N - obj.infCount; %updating susceptible count
            obj.t=obj.t+obj.dt;              %updating current time
        end
    
        function updatedPosition = normalUpdate(obj, curPos, curVel)
            updatedPosition = curPos + obj.dt * curVel;
        end
        
        function updatedVelocity = velocityUpdate(obj, curVel, Force)
            if obj.ACTIVATE_DISTANCING == 1
                reduceMovement = 1.3;
            else
                reduceMovement = 1;
            end
            updatedVelocity = (curVel + obj.dt*Force)/reduceMovement;
        end
        
        function checkBoundaryConditions(obj, posX, posY, index)
            if posX<0   %Checking boundary conditions
                obj.pos(index,1)=obj.pos(index,1)+obj.L;
            end
                
            if posY<0   %Checking boundary conditions
                obj.pos(index,2)=obj.pos(index,2)+obj.L;
            end

            if posX>obj.L   %Checking boundary conditions
                obj.pos(index,1)=obj.pos(index,1)-obj.L;
            end

            if posY>obj.L   %Checking boundary conditions
                obj.pos(index,2)=obj.pos(index,2)-obj.L;
            end
        end
        
        function draw(obj)
            subplot(2,1,1);
            axis equal
            set(gca,'Color','k','XTick',[],'YTick',[]);
            
            counter=1;
            while obj.t<1 % <- /dt is the number of iterations simulation runs
                subplot(2,1,1);
   
                infCountCur = 0;

                obj.step;
                if counter==1
                    for i=1:obj.N %converting points to rounded rectangles
                        obj.phand(i)=rectangle('Position', [obj.pos(i,1)-0.5 obj.pos(i,2)-0.5 0.3 0.3],'Curvature',[1 1],'FaceColor','b','EdgeColor','k');  
                    end
                else
                    for i=1:obj.N %changing position of rectangles on every iteration
                        set(obj.phand(i),'Position', [obj.pos(i,1)-0.5 obj.pos(i,2)-0.5 0.3 0.3]);
                      
                        if (obj.logicalInfected(i) == 1) %Turn dot red if infected
                            set(obj.phand(i),'FaceColor','r');
                        end
                        
                        infStatus = get(obj.phand(i),'FaceColor');
                        if infStatus == [1 0 0] %Call to get FaceColor returns an RGB matrix, so
                                                %red is [1 0 0]
                            infCountCur = infCountCur + 1; %Infection count updates are handled here
                        end
                    end
                end
                obj.infCount = infCountCur;
                counter=counter+1;
                xlim([0 obj.L]);    %setting dimensions of display
                ylim([0 obj.L]);    %setting dimensions of display
                    
                subplot(2,1,2)
                
                plot(obj.timeser, obj.Infected,'r');
                %pause(0.05) %uncomment to make runtime longer, or add to
                %max time also
                drawnow %THIS IS IMPORTANT IT CONTROLS ALL CALLBACKS DONT TOUCH
                
            end
        end
    end
end
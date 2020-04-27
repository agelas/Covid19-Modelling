% this function acts as a regression tool which yields the best fitting
% pair of R0 and recv coefficient values to a data set.
% this is done by running every combination of coefficient values in the
% ODE solving function sirModel, and taking the y values yielded to perform
% a sum of squares test, and determine the best combination by comparing r
% values


% the function takes the a vector of R0s and recvs, the starting y0 values,
% the actual value vector to compare to, and a given interval in the actual cases vector 
% to run the comparison on.

% the function outputs the best fitting coefficients, the R value
% associated with that combination's generated model, and the y values of
% that specific model

function [minR, bestR0, bestRecv, besty] = fitModel(R0s, recvs, y0, startday, finday, cases)
    %set starting values, placeholders to initialize the values
    minR = 0;
    bestR0 = 0;
    bestRecv = 0;
    % find the length of the time interval
    days = finday-startday;
    % initialize array to zeros for storing new y values
    besty = zeros(days, 3);
    for i = 1:max(size(R0s))
        for j = 1:max(size(recvs))
            %iterate thru every combination of coefficients and run it
            %through the sirModel function
            y = sirModel(y0,days+1,R0s(i),recvs(j));
            %find the r value (sum of square of differences) for each
            %iteration
            r = sum((cases(startday:finday)-y(:,2)).^2);
            %check if the generated r value is the smallest generated so
            %far or if it is the first,
            if (minR == 0) || (r < minR)
                %if so, store its associated coefficients and the model y
                %values as the best values so far.
                
                bestR0 = R0s(i);
                bestRecv = recvs(j);
                minR = r;
                besty = y;
            end
        end
    end
    %the function will yield the best values at the end of all
    %iterations
end
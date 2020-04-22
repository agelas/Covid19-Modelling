function [minR, bestR0, bestRecv, besty] = fitModel(R0s, recvs, y0, startday, finday, cases)
    minR = 0;
    bestR0 = 0;
    bestRecv = 0;
    days = finday-startday;
    besty = zeros(days, 1);
    for i = 1:max(size(R0s))
        for j = 1:max(size(recvs))
            y = sirModel(y0,days+1,R0s(i),recvs(j));
            r = sum((cases(startday:finday)-y(:,2)).^2);
            if (minR == 0) || (r < minR)
                bestR0 = R0s(i);
                bestRecv = recvs(j);
                minR = r;
                besty = y;
            end
        end
    end
end
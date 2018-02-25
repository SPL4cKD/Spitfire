function [ tt, tn, solv2] = SatBurn( time, btime, iter)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here


if ~isempty(btime)
tn = find(btime<time(end),1,'last')+1;
tval = [time(1), btime(btime<time(end)), time(end)];
solv2 = 2;
else
    tn = 1;
    tval = [time(1), time(end)];
    solv2 = 1;
end
    for i = 1:tn
        tt{i} = linspace(tval(i),tval(i+1),round((tval(i+1)-tval(i))/(time(end) - time(1)))*iter);
    end
end


    function [k] = TrailSpacing(i, k2, dotNS)
        if i > dotNS(2)*dotNS(1)
            k = i-dotNS(2)*dotNS(1);
        else
            k = k2(i);
        end
    end
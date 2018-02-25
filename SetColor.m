    function [color] = SetColor(BodyN)
        color = [
            249, 219, 026; % sun      1
            122, 120, 122; % mercury  2
            175, 107, 031; % venus    3
            079, 082, 115; % earth    4
            172, 159, 158; % moon     5
            140, 083, 063; % mars     6
            179, 166, 151; % jupiter  7
            206, 172, 117; % saturn   8
            185, 222, 226; % uranus   9
            059, 089, 214; % neptune  10
            171, 135, 112; % pluto    11
            249, 219, 026; % Sat      12
            185, 222, 226; % Sat2     13
            ]/255;
        if BodyN > size(color,1)
            r = randi([0 255],BodyN-size(color,1),3)/255;
            color = vertcat(color, r);
        end
    end
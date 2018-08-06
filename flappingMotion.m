function [L,r,theta] = flappingMotion(L0,Lf,r,theta,steps,i)
    %the cable length and position to generate a flapping motion
    %split the interval into 4 segments
    %1: L0 -> Lf, r = r1
    %2: Lf -> L0, r = r1
    %3: L0 -> Lf, r = r2
    %4: Lf -> L0, r = r2
    theta = theta;
    if i < steps/4
        r = r;
        L = i*(Lf-L0)/(steps/4) + L0;
    elseif i < 2*steps/4
        r = r;
        L = (i-steps/4)*(L0-Lf)/(steps/4) + Lf;
    elseif i < 3*steps/4
        r = -r;
        L = (i-2*steps/4)*(Lf-L0)/(steps/4) + L0;
    else
        r = -r;
        L = (i-3*steps/4)*(L0-Lf)/(steps/4) + Lf;
    end
end
function out = R(x)
    %out = expm(skew(x));
    out = Rx(x(1))*Ry(x(2))*Rz(x(3));
end
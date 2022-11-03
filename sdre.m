function u = sdre(x, sysName, Q, R, A, B)
    if nargin == 4
        A = sys.get_A(x, sysName);
        B = sys.get_B(x, sysName);
    end
    P = icare(A, B, Q, R);
    u = -inv(R) * B' * P * x;
end

function X = accrejrnd(f,g,grnd,c,m,n)
    X = zeros(m,n); % Preallocate memory
    for i = 1:m*n
        accept = false;
        while accept == false
            u = rand();
            v = grnd();
            if c*u <= f(v)/g(v)
               X(i) = v;
               accept = true;
            end
        end
    end
end
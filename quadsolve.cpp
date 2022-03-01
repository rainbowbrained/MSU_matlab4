    X = a: step: b;
    Y = exp(-X.*X); % или другая твоя функция
    Y_primitive = trapz(X, Y); % точное вычисление матлабовской встроенной функцией (вычисляет методом трапеций)
    Y_prim_simp = simpson (X, Y); % метод симпсона
    Y_prim_rect = rectangles (X, Y); % метод прямоугольников
Y_primitive = trapez(X, Y); % Метод трапеций


% task 12: rectangles and simpson methods
% X = [x1, x2, x3, x4, x5] => Y_mid = [f(x2), f(x4)]
% X = [x1, x2, x3, x4] => Y_mid = [f(x2)]
function res = rectangles (X, Y)
    if (length(X) ~= length(Y))
        disp ('Incorrect size grid')
    else 
        n = length(X);
        step = (X(n) - X(1))/n;
        Y_mid = Y([2:2:n]);
        res = 2*step*sum(Y_mid)
    end
end

function res = simpson (X, Y)
    if (length(X) ~= length(Y))
        disp ('Incorrect size grid')
    else 
        n = length(X);
        step = (X(n) - X(1))/(n-1);
        Y_edges = Y([1:2:n]); % a = x1 < x3 < ... < x2n-1 = b
        Y_mid = Y([2:2:n-1]); % x2 < x4 < ... < x2n-2 = b
        res = step*(4*sum(Y_mid) + 2*sum(Y_edges) - Y(n) - Y(1))/3
        if (mod(n, 2) == 0)
            res = res + (Y(n) + Y(n-1))*step/2 
        end
    end
end

function res = trapez (X, Y)
    if (length(X) ~= length(Y))
        disp ('Incorrect size grid')
    else 
        Y_mid = (Y(1:end-1) + Y(2:end))./2;
        n = length(X);
        step = (X(n) - X(1))/n;
        res = step*sum(Y_mid)
    end
end

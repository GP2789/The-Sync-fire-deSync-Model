function [ xy_d ] = phase_compare( x, y_p )
    y_p = cos(y_p) + 1i*sin(y_p); x = conj(cos(x) + 1i*sin(x));
    L = length(y_p); xy_d = zeros(size(x,1), length(L/2 : size(x,2) - L/2));
    for i = L/2 : length(x) - L/2
        %x_p = x(:, i - L/2+1 : i + L/2);
        xy_d(:, i - L/2 + 1) = abs(mean(y_p .* x(:, i - L/2+1 : i + L/2),2));
        %xy_d(i - L/2 + 1) = 1 - abs(sum(exp(1i*(y_p - x_p)))/L);
    end
end


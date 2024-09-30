function alphas = newton_interpolation(x, y)
    if length(x) ~= length(y)
        error('Input vectors must have the same length.');
    end
    n = length(x);
    dd_table = zeros(n, n);
    dd_table(:, 1) = y';
    for j = 2:n
        for i = j:n
            dd_table(i, j) = (dd_table(i, j-1) - dd_table(i-1, j-1)) / (x(i) - x(i-j+1));
        end
    end
    alphas = dd_table(:, 1);
end

function [coefficients] = cubic_spline(x, y)
    n = length(x);
    h = diff(x);
    slope = diff(y) ./ h;
    main_diag = [h(1:end-1) .* 2; 1]; 
    sub_diag = h;
    super_diag = h;
    A = spdiags(sub_diag, main_diag, super_diag); %, -1:1, n-1, n-1)
    b = 6 * diff(slope);
    %[Lcomp, Ucomp, Pcomp, Qcomp] = completePivotingLU(A);
    [Lcomp, Ucomp, Pcomp, Qcomp] = completePivoting2(A);
        %Forward substitution for our given matrix
fs=forwardSubstitution(Lcomp,y);
%Backward substitution for our given matrix
inte_vals=backwardSubstitution(Ucomp,x);
    % Evaluate the interpolating polynomial at x_interpolate
   second_derivatieves=fs;
    
    for i = 1:n-1
        coefficients(i, 1) = (second_derivatives(i+1) - second_derivatives(i)) / (6 * h(i));
        coefficients(i, 2) = second_derivatives(i) / 2;
        coefficients(i, 3) = slope(i) - (h(i) / 6) * (2 * second_derivatives(i) + second_derivatives(i+1));
        coefficients(i, 4) = y(i);
    end
end

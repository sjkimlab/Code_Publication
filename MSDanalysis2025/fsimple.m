function [best_params, fitted_f] = fsimple(x, f_values, D, alpha)
    x = x';
    f_values = f_values'; 
    ft = fittype('log(2*a*exp(x.*b)+c)', 'independent','x', 'coefficients',{'a','b','c'});
    options = fitoptions(ft);
    options.StartPoint = [D, alpha, 0];
    options.Lower = [0, 0, -inf];  % If there is an error change this to [1e-6, 0, -1e-6] to prevent log from going complex.
    options.Upper = [inf, 1, inf];
    best_params = fit(x, f_values, ft, options);
    fitted_f = log(2*best_params.a .* exp(x .* best_params.b) + best_params.c); 
end
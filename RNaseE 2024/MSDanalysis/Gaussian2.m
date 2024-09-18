function y = Gaussian2(x, c, a, mu, sigma, mu2, sigma2)
    %local function for two-population Gaussian
    %c is the scaling factor, mu and mu2 are the means of the first and
    %second Gaussian, sigma and sigma2 are the standard deviations of the
    %first and second Gaussian, a is the fraction of the first Gaussian
    %population, 1-a is the fraction of the second Gaussian population
    y = c*(a/(sigma*sqrt(2*pi))*exp(-0.5*((x-mu)/(sigma)).^2) + (1-a)/(sigma2*sqrt(2*pi))*exp(-0.5*((x-mu2)/(sigma2)).^2));
end
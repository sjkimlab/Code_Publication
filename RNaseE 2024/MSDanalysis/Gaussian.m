function y = Gaussian(x, c, mu, sigma)
    %local function for single-population Gaussian
    %c is the scaling factor, mu is the mean of the Gaussian, sigma is the
    %standard deviation of the Gaussian
   y = c/(sigma*sqrt(2*pi))*exp(-0.5*((x-mu)/(sigma)).^2);
end
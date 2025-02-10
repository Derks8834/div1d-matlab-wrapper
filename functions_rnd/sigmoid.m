function y = sigmoid(x, mu, sigma)
    y = 1./(exp(-1/sigma*(x-mu))+1);
end

function err = errorCoef_vector(z,c)

err = sum( abs(z-c) ) / (size(z,1)*size(c,2));
% err = max( abs(z-c) );
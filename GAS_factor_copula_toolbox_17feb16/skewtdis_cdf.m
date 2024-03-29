function cdf = skewtdis_cdf(x, nu, lambda)
% PURPOSE: returns the cdf at x of Hansen's (1994) 'skewed t' distribution
%---------------------------------------------------
% USAGE: cdf = skewtdis_cdf(x,nu,lambda)
% where: x  = a matrix, vector or scalar 
%        nu = a matrix or scalar degrees of freedom parameter 
%			  lambda = a maxtrix or scalar skewness parameter 
%---------------------------------------------------
% RETURNS:
%        a matrix of cdf at each element of x      
% --------------------------------------------------
% SEE ALSO: tdis_cdf, tdis_rnd, tdis_inv, tdis_prb, skewtdis_pdf
%---------------------------------------------------
%
% Based on tdis_cdf.m from the "Spatial Econometrics"
% toolbox of James P. LeSage
% http://www.spatial-econometrics.com/
%
%  Andrew Patton
%
%  25 June, 2001

% This code was used in: 
%
% Patton, Andrew J., 2004, "On the Out-of-Sample 
% Importance of Skewness and Asymmetric Dependence for
% Asset Allocation", Journal of Financial Econometrics, 2(1), 130-168.


[T,k] = size(x);
if size(nu,1)<T;
   nu = nu(1)*ones(T,1);
end
if size(lambda,1)<T;
   lambda = lambda(1)*ones(T,1);
end
c = gamma((nu+1)/2)./(sqrt(pi*(nu-2)).*gamma(nu/2));
a = 4*lambda.*c.*((nu-2)./(nu-1));
b = sqrt(1 + 3*lambda.^2 - a.^2);

y1 = (b.*x+a)./(1-lambda).*sqrt(nu./(nu-2));		
y2 = (b.*x+a)./(1+lambda).*sqrt(nu./(nu-2));

 cdf = (1-lambda).*tcdf(y1,nu).*(x<-a./b);
 cdf = cdf + (x>=-a./b).*((1-lambda)/2 + (1+lambda).*(tcdf(y2,nu)-0.5));

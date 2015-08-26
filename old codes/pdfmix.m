function [out] = pdfmix(x,mu1,mu2,Sigma1,Sigma2,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
z=(Sigma1)\(x-mu1)';
w=(Sigma2)\(x-mu2)';
d=length(x);
mu=zeros(1,d);
Sigma=eye(d);
out=p*mvnpdf(z,mu,Sigma)+ (1-p)*mvnpdf(w,mu,Sigma);
end

%addpath('C:\Documents and Settings\Admin\My Documents\MATLAB\MC3');

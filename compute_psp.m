function [ EPSP ] = compute_psp( tau_syn, T )
%COMPUTE_PSP Summary of this function goes here
%   Detailed explanation goes here
EPSP = (exp(1)*T./tau_syn).*exp(-T./tau_syn).*heaviside(T); 

end


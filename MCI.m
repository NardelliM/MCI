function [Lambda,mci_sh,mci_lo]=MCI(X,m,r,n,tau,scl,varargin)
%
% MCI: Multichannel Complexity Index
% This function implements the MCI algorithm described in "Multichannel 
% Complexity Index (MCI) for a Multi-Organ Physiological Complexity Assessment", 
% Physica A, 2019
% MCI relies on a novel method for the reconstruction of the multivariate 
% phase space, where each series is embedded using its proper time delay. MCI
% accounts for the estimation of phase space distances using fuzzy rules, 
% and may be computed at two different ranges of time-scale values to 
% investigate short- and long-term dynamics.
%
% Syntax:
%   [Lambda,mci_sh,mci_lo]=MCI(X,m,r,n,tau,scl,scales_sh,scales_lo]
%
%   Inputs:
% X: matrix related to the c-variate time series - a matrix of size c (the 
% number of channels) x N (the number of sample points for each channel)
% m: vector of embedding dimensions - size 1 x c, following the order of channels
% r: scalar value related to the threshold value (usually equal to 0.15)
% n: fuzzy power (usually equal to 2)
% tau: vector of time delays - size 1 x c, following the order of channels
% scl: number of scales to be analysed 
% scales_sh: vectors of short-time scales for the computation of mci_sh (default [1:5])
% scales_lo: vectors of long-time scales for the computation of mci_lo (default [6:scl])
%
%   Output:
% Lambda: values of multivariate fuzzy entropy for each scale factor (vector 
% of size 1 x scales)
% mci_sh: multichannel complexity index computed as the area under the curve 
% of Lambda for short scales (default from 1 to 5). If scl<=5, mci_sh is computed 
% considering the scales lower than scl.
% mci_lo: multichannel complexity index computed as the area under the curve 
% of Lambda for long scales (default from 6 to the highest). If scl<=5, mci_lo=[].
%
% Ref:
% [1] M. Nardelli, E.P.Scilingo, G.Valenza "Multichannel Complexity Index (MCI) 
%     for a Multi-Organ Physiological Complexity Assessment", Physica A: Statistical
%     Mechanics and its Applications, 2019, https://doi.org/10.1016/j.physa.2019.121543
%
% If you use the code, please make sure that you cite reference [1].
%
% You may contact the author by e-mail: 
% Mimma Nardelli
% m.nardelli@ing.unipi.it
%
% Copyright Mimma Nardelli, Enzo Pasquale Scilingo, and Gaetano Valenza,
% May 2019
%______________________________________________________________________________
%
% File:                         MCI.m
% Last revised:                 28 May 2019
% ______________________________________________________________________________
%
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the Creative Commons Attribution 4.0 for more details.
%
% ______________________________________________________________________________

if nargin<6
    error('Incorrect call to MCI: missing inputs');
end

if scl>5
    argin = {1:5, 6:scl};
else
    argin = {1:scl};
end
i = ~cellfun(@isempty, varargin);
argin(i) = varargin(i);
[scales_sh, scales_lo] = deal(argin{:});

if scl<scales_sh(end)
    error('Incorrect call to MCI: scl is lower than scales_sh upper limit');
end

if scl<scales_lo(end)
    error('Incorrect call to MCI: scl is lower than scales_lo upper limit');
end

% Multichannel signals may have different amplitude ranges: a normalization 
% to unit standard deviation is computed using zscore

X=zscore(X');
X=X';
N=size(X,2);C=size(X,1);

% Coarse-graining procedure: for each value of scale factor, the matrix of
% scaled series Xsc is computed. Each scaled series is obtained by averaging 
% the data points within non-overlapping windows

for s=1:scl
    Xsc=zeros(C,floor(size(X,2)/s));
    len_c=zeros(1,C);
    for c=1:size(X,1)
        len_c(c)=length(1:tau(c):floor(size(X,2)/s));
        Xsc(c,1)=mean(X(c,1:s));
        for i=1:((size(X,2)/s)-1)
            Xsc(c,i+1)=mean(X(c,1+s*i:s*i+s));
        end
    end
    
% Multivariate phase space reconstruction, following the MCI estimation procedure [1]

    l_min= min(len_c);
    Xsc_tau=zeros(c,l_min);
    
    for c=1:size(X,1)
        Xtau=Xsc(c,1:tau(c):end);
        Xsc_tau(c,:)=Xtau(1:l_min);
    end
    Xsc_median=median(Xsc_tau);
    
% Fuzzy Entropy computation     
    
    Lambda(s)=fuzzyen(Xsc_median,max(m),r*sum(std(Xsc')),n,1);
    
end

% MCI values are calculated as the areas under the curve of Lambda values
% as a function of the scale factors, for short and long scales 

mci_sh=trapz(Lambda(scales_sh));
if scl>scales_sh(end)
    mci_lo=trapz(Lambda(scales_lo));
else
    mci_lo=[];
end


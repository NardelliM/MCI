function [fuz_en]=fuzzyen(X,m,rf,fp,tau)
% This function compute Fuzzy Entropy of the time series in vector X.
% The computation is based on the use of an the exponential membership
% function as in:
% [1] M. Nardelli, E.P.Scilingo, G.Valenza "Multichannel Complexity Index 
%    (MCI) for a Multi-Organ Physiological Complexity Assessment", Physica A, 2019.
%     https://doi.org/10.1016/j.physa.2019.121543
% [2] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized 
%     Multiscale Fuzzy Entropy: A Tool for Complexity Analysis of Multichannel 
%     Signals", Physica A, 2016.
% [3] W. Chen, Z. Wang, H. Xie, and W. Yu,"Characterization of surface EMG 
%     signal based on fuzzy entropy", IEEE Transactions on neural systems 
%     and rehabilitation engineering, vol. 15, no. 2, pp.266-272, 2007.
%
%
% If you use the code, please make sure that you cite references [1],[2],[3].
%
% m.nardelli@ing.unipi.it
%
% Input:
% X: time series vector of size 1xN, where N is the number of samples
% m: embedding dimension
% rf:radius used for the comparison of embedding vectors (usally rf=0.15*std(X))
% fp=fuzzy power (usually fp=2)
% tau=time delay
%
% Output
% fuz_en:fuzzy entropy related to time series X

% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the Creative Commons Attribution 4.0 for more details.


%% Starting from the embedding dimension m:
n=m*tau;
N=size(X,2)-n;

% All the embedded vectors are created based on the Takens embedding
% Theorem, using the input values related to m and tau

vett1=[];
parfor i=1:N
t1(i,:)=X(i:tau:i+(m-1)*tau);
end
vett1=horzcat(vett1,t1);
dis=pdist(vett1,'chebychev'); % Chebychev distance is computed between all 
%possible pairs of embedding vectors with dimension m
gamma=exp((-dis.^fp)/rf); % The the exponential membership function gamma is computed
fuz_m1=sum(gamma)*2/(N*(N-1));


%% Increasing the embedding dimension to m+1
m=m+1;
n=m*tau;
N=size(X,2)-n;
vett2=[];
parfor i=1:N
t2(i,:)=X(i:tau:i+(m-1)*tau);
end
vett2=horzcat(vett2,t2);
dis2=pdist(vett2,'chebychev'); %Chebychev distance is computed between all 
% possible pairs of embedding vectors with dimension m+1

gamma2=exp((-dis2.^fp)/rf); % The the exponential membership function gamma is computed

fuz_m2=sum(gamma2)*2/(N*(N-1));

%% calculating Fuzzy Entropy

fuz_en=log(fuz_m1/fuz_m2);

end


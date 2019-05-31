# MCI
Algorithm for the computation of Multichannel Complexity Index (MCI)

This code implements the MCI algorithm for the computation of complexity 
index in multichannel system, described in:

[1] M. Nardelli, E.P. Scilingo, G. Valenza "Multichannel Complexity Index (MCI)
for a Multi-Organ Physiological Complexity Assessment", Physica A: Statistical
Mechanics and its Applications, 2019, https://doi.org/10.1016/j.physa.2019.121543

This metrics relies on a novel method for the reconstruction of the
multivariate phase space, where each series is embedded using its proper 
time delay. Then, MCI accounts for the estimation of phase space distances 
using fuzzy rules, and may be computed at two different ranges of time-scale
values to investigate short- and long-term dynamics.

The algorithm was evaluated using three-channel white gaussian noise and 1/f 
noise systems, synthetic series from the Henon map and Rossler attractor, 
and publicly-available physiological series (see the article). 

___________________________________________________________________________

Copyright (C) 2019 Mimma Nardelli, Enzo Pasquale Scilingo, Gaetano Valenza

This program is free software; you can use it under the terms of the 
Creative Commons License: Attribution 4.0 International. 

If you use this program in support of published research, please include a
citation of the reference [1]. If you use this code in a software package,
please explicitly inform the end users of this copyright notice and ask them
to cite the reference above in their published research.
___________________________________________________________________________

To use this software from Matlab, simply call the MCI function in the 
path/folder. Type 'help MCI' from the Matlab command window for help on the 
command's syntax and input/output arguments.

The software does not come with a GUI. Assuming that:
-'X' is a matrix related to the c-variate time series -a matrix of size c 
(the number of channels) x N (the number of sample points for each channel),

-'m' is the vector of embedding dimensions,

-'r' is the scalar value related to the threshold value (it is usually equal
to 0.15),

-'n' is the fuzzy power (it is usually equal to 2),

-'tau' is the vector of time delays,

-'scl' is the number of scales studied (usually 20),


the following example performs the MCI analysis (with default parameters) 
and plots the results:

[Lambda,mci_sh,mci_lo]=MCI(X,m,0.15,2,tau,20);

figure,

plot([1:20],Lambda,'LineWidth',1.5,'Color','k');

ylabel('\Lambda _{MCI}','Fontsize',14)

xlabel('Scale factor (\beta)','Fontsize',14)

str = strcat('MCI_{short}= ',num2str(mci_sh));

str2 = strcat('MCI_{long}= ',num2str(mci_lo));

ylim=get(gca,'ytick');

xlim=get(gca,'xtick');

text(xlim(1),ylim(end-1),{str,str2},'BackgroundColor',[0.98 0.98 0.98],'Fontsize',13)





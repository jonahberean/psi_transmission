# psi_transmission

## plan of work 
1. By going through the eLog, and documentation in the .txt files, a map of the different set-ups and storage times should be created. We want the different files to be renamed according to some convention that allows analysis over many files belonging to the same configurations.
2. The time histograms for individual runs, i.e. individual .tof files, show an initial peak due to target irradiation, we need to appropriately select an integration window that excludes this initial peak.
3. Produce analysis that shows neutron counts vs. storage time. This will allow the determination of the pre-storage lifetime (a characteristic decay time) from an exponential fit.
4. These results can be compared to Stewart's prestorage simulations, so that the imaginary Fermi potential may be determined.
5. Does this data include storage lifetime measurements of the guides that were used? If so, we can then use simulations to determine the specular v. diffuse reflection probability.
6. Using the runs that act as a sort of 'normalization monitor', the performance of the source can be fit to some decay. Then the above analysis can be renormalized using this decay.
7. Finally, determine transmission (a more nuanced process) as a function of storage time.

## some organizational tasks
1. arrange misc/ better

## notes

There are three different data sets with proton beam data. (!!! which one was used, where are they all from?)

This project has been set up using PyScaffold 3.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.

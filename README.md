## ASAF
Codes for Stochastic GW Background Analysis for the project done with Prof. Sanjit Mitra at IUCAA.
The module contains the functions I created to find order 1, 2, 3, and 4 peaks in ASAF maps.
O1 peaks are found by comparing the SNR of a pixel to its nearest neighbours. If its SNR is greater than all its neighbours, it is counted as a peak. O2 peaks are found by comparing a pixel to its nearest and next to nearest neighbours; and so on for O3 and O4 peaks...

The function get_peaks_new finds the peaks for a map. A list of neighbours need to passed onto this function as input. This list is found by the functions like get_o1neighbours, get_o2neighbours and so on, depending on the order of the peaks we want. These functions return the healpy indices of the all the neighbours of the all the pixels in a map of given $N_{side}$.

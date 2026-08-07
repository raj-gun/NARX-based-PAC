# Filtering-based Phase-Amplitude Coupling Methods

This folder contains MATLAB implementations of the filtering-based phase-amplitude coupling (PAC) methods used for comparison with the NARX-based PAC method.

## Implemented Methods

The following PAC measures are included:

- **Modulation Index (MI)** proposed by Tort et al. [1].
- **Mean Vector Length (MVL)** proposed by Canolty et al. [2].
- **Direct PAC estimator** proposed by Özkurt and Schnitzler [3].
- **General Linear Model (GLM) PAC estimator** proposed by Penny et al. [4].

These methods are used for the comparisons presented in the paper and supplementary figures.

## Original Source Code

The implementations in this folder were adapted from:

- MATLAB code accompanying Özkurt and Schnitzler [3], which was adapted from code written for Kramer et al.
- The [TortLab phase-amplitude-coupling repository](https://github.com/tortlab/phase-amplitude-coupling), containing MATLAB routines for computing the Modulation Index and comodulograms described by Tort et al. [1].

For the comparisons used in the NARX-based PAC study, the filtering bandwidths were made adaptive following Berman et al. [5].

## References

[1] A. B. L. Tort, R. Komorowski, H. Eichenbaum, and N. Kopell, “Measuring phase-amplitude coupling between neuronal oscillations of different frequencies,” *Journal of Neurophysiology*, vol. 104, no. 2, pp. 1195–1210, 2010. doi:10.1152/jn.00106.2010.

[2] R. T. Canolty, E. Edwards, S. S. Dalal, M. Soltani, S. S. Nagarajan, H. E. Kirsch, M. S. Berger, N. M. Barbaro, and R. T. Knight, “High gamma power is phase-locked to theta oscillations in human neocortex,” *Science*, vol. 313, no. 5793, pp. 1626–1628, 2006. doi:10.1126/science.1128115.

[3] T. E. Özkurt and A. Schnitzler, “A critical note on the definition of phase–amplitude cross-frequency coupling,” *Journal of Neuroscience Methods*, vol. 201, no. 2, pp. 438–443, 2011. doi:10.1016/j.jneumeth.2011.08.014.

[4] W. D. Penny, E. Duzel, K. J. Miller, and J. G. Ojemann, “Testing for nested oscillation,” *Journal of Neuroscience Methods*, vol. 174, no. 1, pp. 50–61, 2008. doi:10.1016/j.jneumeth.2008.06.035.

[5] J. I. Berman, J. McDaniel, S. Liu, L. Cornew, W. Gaetz, T. P. L. Roberts, and J. C. Edgar, “Variable bandwidth filtering for improved sensitivity of cross-frequency coupling metrics,” *Brain Connectivity*, vol. 2, no. 3, pp. 155–163, 2012. doi:10.1089/brain.2012.0085.

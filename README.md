# Grassmannian Constellations for MIMO Noncoherent Communications

Structure family of Grassmannian constellations for MIMO noncoherent communications. The proposed constellation design is built upon the geodesic curves of the Grassmann manifold, thereby exploiting its underlying geometric structure. Remarkably, the obtained space-time constellations have a single nonzero entry per row, meaning that a single antenna is active in transmission per time slot.

A preprint version of the paper explaining this design can be found at: https://arxiv.org/abs/2510.15070

Contributors:
- [Álvaro Pendás-Recondo](https://scholar.google.com/citations?user=akuJpAIAAAAJ&hl=en)
- [Enrique Pendás-Recondo](https://scholar.google.com/citations?user=hlbBliwAAAAJ&hl=en)


This repository contains examples of constellations obtained with this approach and the code to generate them. Variables are:

- T -> Duration of the coherence block in symbols.
- M -> Number of transmit antennas.
- L -> Constellation size or number of points.
- N -> Number of receive antennas.

This code is intended for research purposes. Some tasks could be further automated; for example, the selection of geodesic indices is currently performed by running the script diametral_set.m and then manually adjusting the indices.

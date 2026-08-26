# Frontal_Gradients_Boundaries

This repository contains code used in the following preprint:

   [not yet submitted to bioRciv]

This repository contains code for preprocessing, analysing, and making figures for the frontal gradients boundaries.

Running most of this code requires access to the data on the SCC.

Organization:

    scripts: contains key scripts for file i/o and organization, file conversion, preprocessing, ROI creation, and analysis
          >GLM_scripts: scripts creating some of the extra GLM contrasts used in analysis
          >ROI_creation: script for creating the individual subj-level ROIs used in analysis
          >ROI_probabilistics_scripts: scripts for permuting subj-level ROI labels and combining into probabilistics used in permutation testing for COM shift analysis. 
          >analysis: all analyses for this project (COM shift, domain general dice, gradient/boundary modeling, behavioral accuracy, etc)
          >misc: file reorg scrip and screenshotting scripts for freeview
          >plotting: scripts that only plot data used in figures


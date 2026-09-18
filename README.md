# Astrocyte–GABAergic Neuron Metabolic Interaction Model

This repository contains the computational models, scripts, and supporting files used to investigate metabolic alterations in astrocytes and GABAergic neurons across the cognitive impairment continuum, from control conditions to mild cognitive impairment (MCI) and Alzheimer's disease (AD).

The study integrates transcriptomic information with genome-scale metabolic models (GEMs) to characterize cell-specific metabolic states and the metabolic interactions between astrocytes and GABAergic neurons.

## Overview

Astrocytes and neurons maintain a close metabolic relationship that is essential for neuronal function and brain homeostasis. Alterations in this metabolic coupling may contribute to the metabolic vulnerability observed during cognitive decline and Alzheimer's disease progression.

The computational workflow implemented in this repository includes:

1. Integration of transcriptomic information into cell-specific genome-scale metabolic models.
2. Reconstruction and analysis of astrocyte and GABAergic neuron metabolic models.
3. Constraint-based characterization of metabolic states across disease conditions.
4. Flux Balance Analysis (FBA) and Flux Variability Analysis (FVA).
5. Construction of an astrocyte–GABAergic neuron community metabolic model.
6. Analysis of metabolite exchange and metabolic interactions between both cell types.
7. Comparison of metabolic behavior across the Control–MCI–AD continuum.

## Study Conditions

The metabolic models represent four stages of cognitive impairment:

- Control
- Early mild cognitive impairment (E-MCI)
- Advanced mild cognitive impairment (A-MCI)
- Alzheimer's disease (AD)

These conditions were analyzed independently at the cell-specific level and subsequently integrated into astrocyte–neuron community models.

## Computational Workflow

The general workflow of the study is:

Human hippocampal transcriptomic data
                │
                ▼
Transcriptomic deconvolution
                │
                ▼
Cell-specific expression profiles
        ┌───────┴────────┐
        ▼                ▼
    Astrocytes      GABAergic neurons
        │                │
        ▼                ▼
 Cell-specific genome-scale
      metabolic models
        │                │
        └───────┬────────┘
                ▼
      Constraint-based analysis
           (FBA / FVA)
                │
                ▼
 Astrocyte–GABAergic neuron
     community metabolic model
                │
                ▼
 Metabolic exchange and interaction
                │
                ▼
 Comparison across
 Control → E-MCI → A-MCI → AD

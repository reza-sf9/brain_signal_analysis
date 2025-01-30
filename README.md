# Brain Signal Analysis

This repository contains three different projects related to brain signal analysis using coherence modeling and switching dynamics.

## 1. Latent Dynamical Coherence Model

### Overview
In our previous work, we developed a latent dynamical modeling framework called **state-space global coherence** to capture slow-changing dynamics in network-level coherence. In this research, we extend this approach by developing a more general class of **state-space coherence models (SS-Coh)** that can capture **fast and switching changes** in network-level rhythmic dynamics.

### Data
We use **EEG data from human patients under general anesthesia** collected in **Emery Brown’s laboratory**. The dataset records 64-channel EEG at a 250 Hz sampling rate during **induction and emergence from propofol anesthesia**.

### Results
LDCM incorporates **switching mechanisms** to improve estimation accuracy over SS-Coh models.

![LDCM Results](path/to/ldcm_results.png) <!-- Placeholder for LDCM results image -->

## 2. State Space Coherence (SS-Coh)

### Overview
State Space Coherence (**SS-Coh**) provides a **robust algorithm** for estimating **Global Coherence (GCoh)**. This approach improves tracking of functional connectivity changes during different states of consciousness.

### Data
Similar to LDCM, SS-Coh utilizes **EEG data from patients under general anesthesia**. The dataset is publicly available from **MGH Human Research Committee**.

### Results
Applying SS-GCoh to anesthesia data shows **functional circuit changes associated with consciousness states**.

![SS-Coh Results](path/to/sscoh_results.png) <!-- Placeholder for SS-Coh results image -->

### Implementation
To reproduce these results, run the provided **GMM code** to obtain functional modes corresponding to different consciousness states.

## 3. Switching Dynamic

### Overview
This project evaluates **switching mechanisms in brain activity** using **Global Coherence Algorithm**. We investigate how **neural circuits engage differently** across cognitive tasks.

### Data
EEG data was recorded in **14-minute sessions** using a **20-channel dry-electrode EEG**. Subjects performed cognitive tasks including:
- **Eyes Open (EO)**
- **Reading (RD)**
- **Mental Subtraction (MS)**
- **Eyes Closed (EC)**

### Results
Task-switching analysis reveals that distinct **functional circuits emerge during different cognitive tasks**.

![Switching Dynamic Results](path/to/switching_results.png) <!-- Placeholder for switching dynamic results image -->

### Implementation
Run the **main code** to obtain **task-specific functional clusters**.

## Contributors
- **Reza Saadati Fard**
- **Collaborators from Emery Brown’s Laboratory**

## References
1. Purdon et al., "Anesthesia and EEG coherence," *Nature Neuroscience*, 2019.
2. Brown et al., "State-space modeling in brain signals," *Journal of Neuroscience Methods*, 2020.

# 🧬 HIV-1 env Gene Evolutionary Analysis

## 📌 Overview

This project presents a bioinformatics analysis of sequence variation in the HIV-1 **env gene** across geographically distinct strains. The study investigates how nucleotide-level variation reflects evolutionary pressures such as structural constraint and immune evasion.

---

## 🎯 Research Question

How does nucleotide sequence variation in the HIV-1 env gene reflect mechanisms of host immune evasion?

---

## 🧠 Biological Background

The HIV-1 env gene encodes envelope proteins (gp120 and gp41) responsible for:

* Host cell entry (CD4 and co-receptors)
* Immune system recognition

The env gene is known for:

* High variability → immune escape
* Functional constraints → structural conservation

---

## 🧪 Data

* 4 HIV-1 env gene sequences:

  * Nigeria
  * USA
  * Zambia
  * South Africa
* Source: NCBI GenBank / HIV databases
* Data type: DNA sequences (aligned)

---

## ⚙️ Methodology

### 1. Multiple Sequence Alignment

Sequences were aligned using MUSCLE to identify conserved and variable regions.

### 2. Position-wise Analysis

For each nucleotide position, the following metrics were computed:

* **Conservation Score**

  * Frequency of the most common nucleotide
* **Shannon Entropy**

  * Measures sequence diversity
* **Gap Fraction**

  * Frequency of insertions/deletions

### 3. Classification

Positions were categorized using **quantile-based thresholds**:

* Highly conserved
* Moderately conserved
* Variable

---

## 📊 Results

### Key Observations

* **High Conservation (~0.8-1.0)**

  * Indicates strong purifying selection
  * Suggests structural and functional importance

* **Moderate Variability**

  * Regions tolerate mutations
  * Potential adaptive flexibility

* **High Entropy + Gap Regions**

  * Localized mutation hotspots
  * Likely associated with immune evasion

---

## 🧬 Biological Interpretation

The results support a structured evolutionary model:

> The HIV-1 env gene exhibits a conserved functional backbone with localized regions of variability, reflecting a balance between structural constraint and immune-driven evolution.

---

## ⚠️ Limitations

* Analysis based on a **partial env gene fragment**
* Small dataset (**n = 4 sequences**)
* Approximate mapping to protein functional regions
* No phylogenetic or codon-level selection analysis

---

## 💻 Tools & Libraries

* Biopython
* MUSCLE (alignment)
* pandas
* matplotlib
* NumPy

---

## 📁 Outputs

* `env_alignment.aln` → aligned sequences
* `results.csv` → position-wise metrics
* Evolutionary profile plot (conservation, entropy, gap fraction)
* Python pipeline script

---

## 🚀 Key Skills Demonstrated

* Bioinformatics pipeline development
* Sequence alignment and analysis
* Evolutionary interpretation
* Data visualization
* Scientific reasoning

---

## 📌 Conclusion

This project demonstrates how sequence-level variation in viral genes can be quantitatively analyzed to reveal underlying evolutionary mechanisms, providing insight into host-pathogen interactions.

---

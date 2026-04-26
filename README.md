🧬 HIV-1 env Gene Evolutionary Analysis
📌 Overview
This project presents a bioinformatics analysis of sequence variation in the HIV-1 env gene across geographically distinct strains. The study investigates how nucleotide-level variation reflects evolutionary pressures such as structural constraint and immune evasion.

🎯 Research Question
How does nucleotide sequence variation in the HIV-1 env gene reflect mechanisms of host immune evasion?

🧠 Biological Background
The HIV-1 env gene encodes envelope proteins (gp120 and gp41) responsible for:
•	Host cell entry (CD4 and co-receptors)
•	Immune system recognition
The env gene is known for:
•	High variability → immune escape
•	Functional constraints → structural conservation

🧪 Data
•	4 HIV-1 env gene sequences:
o	Nigeria
o	USA
o	Zambia
o	South Africa
•	Source: NCBI GenBank / HIV databases
•	Data type: DNA sequences (aligned)

⚙️ Methodology
1. Multiple Sequence Alignment
Sequences were aligned using MUSCLE to identify conserved and variable regions.
2. Position-wise Analysis
For each nucleotide position, the following metrics were computed:
•	Conservation Score
o	Frequency of the most common nucleotide
•	Shannon Entropy
o	Measures sequence diversity
•	Gap Fraction
o	Frequency of insertions/deletions
3. Classification
Positions were categorized using quantile-based thresholds:
•	Highly conserved
•	Moderately conserved
•	Variable

📊 Results
Key Observations
•	High Conservation (~0.8-1.0)
o	Indicates strong purifying selection
o	Suggests structural and functional importance
•	Moderate Variability
o	Regions tolerate mutations
o	Potential adaptive flexibility
•	High Entropy + Gap Regions
o	Localized mutation hotspots
o	Likely associated with immune evasion

🧬 Biological Interpretation
The results support a structured evolutionary model:
The HIV-1 env gene exhibits a conserved functional backbone with localized regions of variability, reflecting a balance between structural constraint and immune-driven evolution.

⚠️ Limitations
•	Analysis based on a partial env gene fragment
•	Small dataset (n = 4 sequences)
•	Approximate mapping to protein functional regions
•	No phylogenetic or codon-level selection analysis

💻 Tools & Libraries
•	Biopython
•	MUSCLE (alignment)
•	pandas
•	matplotlib
•	NumPy

📁 Outputs
•	env_alignment.aln → aligned sequences
•	results.csv → position-wise metrics
•	Evolutionary profile plot (conservation, entropy, gap fraction)
•	Python pipeline script


🚀 Key Skills Demonstrated
•	Bioinformatics pipeline development
•	Sequence alignment and analysis
•	Evolutionary interpretation
•	Data visualization
•	Scientific reasoning

📌 Conclusion
This project demonstrates how sequence-level variation in viral genes can be quantitatively analyzed to reveal underlying evolutionary mechanisms, providing insight into host-pathogen interactions.



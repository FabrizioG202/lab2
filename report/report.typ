#set page(
  columns: 2,
)

#set text(
  size: 11pt,
)

#set par(
  leading: 0.56em,
  justify: true,
)

#place(
  top + center,
  scope: "parent",
  float: true,
  text(1.4em, weight: "bold")[
    A comparative analysis of Von-Heijne, SVM and Perceptron Models for the detecton of Signal Peptides.
  ],
)

= Abstract
#text(weight: "bold")[Motivation:] This is the motivation

#text(weight: "bold")[Results:] This are the results.

= Introduction
For a protein to enter the secretory pathway, in both eukaryotic and prokaryotic cells, it must be endowed with a specific target signal. Often, this signal takes the shape of a short sequences located at the N-terminus of proteins #link("http://www.ncbi.nlm.nih.gov/pubmed/2197415")[ref]. The signal peptide (SP) has a distinct three-domain structure, depicted in @sp_structure, with a positively charged N-terminal region (n-region), a central hydrophobic region (H-region) and a more polar C-terminal region (C-region) containing the cleavage site #link("http://www.ncbi.nlm.nih.gov/pubmed/2197415")[ref]. Given their importance in many aspects of cell biology, the accurate detection of SPs is a crucial task in bioinformatics, which has been tackled before with a variety of approaches, including machine learning (SVMs #link("https://pubmed.ncbi.nlm.nih.gov/19470175")[ref], #link("https://pubmed.ncbi.nlm.nih.gov/21959131/")[ref] and Bayesian networks #link("http://www.ncbi.nlm.nih.gov/pubmed/18989393")[ref], Hidden Markov Models #link("https://doi.org/10.1016/j.jmb.2004.03.016")[ref]), and, more recently, deep learning #link("https://academic.oup.com/bioinformatics/article/34/10/1690/4769493?login=false")[ref].

#figure(
  image(".imgs/sp_structure.png"),
  caption: [Signal peptide structure],
) <sp_structure>


With this work, we propose a comparative analysis of two different approaches for the detection of SPs, namely the Von-Heijne model a Support Vector Machine (SVM) based method.


= Materials and Methods
== Data Collection
Reviewed non-fragment protein sequences from the eukaryota taxa were queried from UniprotKB, and collected into 2 distinct sets, positive and negative.
For Positive sequences, we selected sequences having experimental evidence for the presence of a signal peptide (`ft_signal_exp:*`)
For the negative set, proteins not annotated with a signal peptide and endowed with experimental evidence for subcellular localization in non-secretory comparments, namely nucleus (`cc_scl_term_exp:SL-0191`), peroxisome (`cc_scl_term_exp:SL-0204`), cell membrane / plasma membrane (`cc_scl_term_exp:SL-0039`), cytosol (`cc_scl_term_exp:SL-0091`), plastid (`cc_scl_term_exp:SL-0209`) or mitochondrion (`cc_scl_term_exp:SL-0173`). For sequences in the negative set, we gathered information about the presence of a transmembrane (TM) helix, whose  hydrophobicity profile mimicks the one of a signal peptide.

Further selection on the positive entries was carried out in order to select proteins which had informaton about the cleavage site and a Signal Peptide with a length of at least 14bp.

== Clustering

Sequences were clustered using mmseqs2 #link("https://www.nature.com/articles/nbt.3988")[ref], using parameters `--min-seq-id 0.3, cluster-mode 1 -cov-mode 0 -c 0.4`, in order to cluster sequences with a similarity of at least 30\% and an alignment spanning at least 40\% of the total sequence length. The number of sequences in the datasets is reported in @seq-counts.

#figure(
  table(
    columns: 3,
    align: (left, center, center),
    [*Dataset*], [*Before Clustering*], [*After Clustering*],
    [Positive], [2942], [1093],
    [Negative (total)], [20806], [9028],
    [Negative (w/o TM)], [18301], [8113],
    [Negative (with TM)], [2505], [915],
  ),
  caption: [Dataset sizes before and after clustering],
)<seq-counts>

== Von-Heijne Model

The full dataset was split into a training and a test set with a ratio of 0.8. A Position-Specific Weight Matrix (PSWM) was then generated using 16bp windows (13bp upstream to 2bp downstream) surrounding the cleavage site. Amino acid frequencies were normalized using background amino acid frequencies from SwissProt. After initializing the matrix with pseudocounts of 1, the weights are computed as follows:

$
  w(a, i) = log_2(f(a, i) / f(a))
$

Where $f(a, i)$ is the frequency of amino acid $a$ at position $i$ in the training set, and $f(a)$ is the background frequency of amino acid $a$ in SwissProt. The score for a given sequence is then computed as:
$
  S = max_(j=1)^(90-16+1) sum_(i=1)^(16) w(a_(j+i-1), i)
$

scores are computed over the first 90bp of the given sequence, using a sliding window 16bp long. Then, a sequence was labeled as positively endowed with a SP if the Score $S$ is greater than a threshold $t$.

$
  "label" = cases(
    "positive" & "if" S > t,
    "negative" & "if" S <= t
  )
$

The optimal threshold for the dataset was computed using 5-fold cross validation.

== Feature Extraction
To build a model that can deal with sequences, they must be encoded into numerical representation. There are many possibility of extracting meaningful numerical features from a sequence, for the SVM and MLP model, we employed 2 kinds of features:
- Aminoacid composition up to a cutoff $k$ (`n_term_composition`) and from $k$ to an hardcoded limit $N$ (set to $90$) (`c_term_composition`)
- Protein Scales computed on the first $N$ aminoacids with a window size of $5$ (@feature-example). An Hydrophobicity scale (cit), alpha-helix tendency scale (cit), transmembrane tendency (cit), bulkiness (cit), and polarity (cit) were used.

For this experiment, $k$ was chosen to be the average position of the cut site, roughly 22 aminoacids from the N-Terminus.

#figure(
  image(".imgs/feature_extraction_example.svg"),
  caption: [Example of feature extraction for a given sequence. Protein scales are computed on the first $N$ aminoacids with a window size of 5],
) <feature-example>


== Support Vector Machine
We used a grid search to find optimal hyperparameters, testing both RBF and linear kernels. For the RBF kernel, we evaluated C values of 0.1, 1, 2, 4, and 8, combined with gamma values of "scale", 0.001, 0.01, 0.1, 1, and 2. For the linear kernel, we tested C values of 0.1, 1, 2, 4, and 8. Each model was tested on a and evaluated on an 80\% training split. The model yielding the highest f1 score was found to have a `rbf` kernel, with $C=2$ and $gamma=0.01$, with a mean f1 score across the 5 folds ran on the training set of $0.87$.


== Evaluation and comparison
To evaluate and compare the performance of the models, we computed several metrics on the test set. The metrics were calculated based on the following formulas:

*Accuracy:* The proportion of correctly classified instances among all instances.
$
  "Accuracy" = (T P + T N) / (T P + T N + F P + F N)
$

*Precision:* The proportion of true positive predictions among all positive predictions.
$
  "Precision" = (T P) / (T P + F P)
$

*Recall (Sensitivity):* The proportion of true positive instances that were correctly identified.
$
  "Recall" = (T P) / (T P + F N)
$

*F1 Score:* The harmonic mean of precision and recall.
$
  #text(size: 8pt)[$F_1 = 2 dot (("Precision" dot "Recall")) / ("Precision" + "Recall") = (2 T P) / (2 T P + F P + F N)$]
$

*Matthews Correlation Coefficient (MCC):* A balanced measure that takes into account all four confusion matrix categories.
$
  #text(size: 8pt)[$"MCC" = (T P dot T N - F P dot F N) / sqrt((T P + F P)(T P + F N)(T N + F P)(T N + F N))$]
$

Where $T P$ (True Positives), $T N$ (True Negatives), $F P$ (False Positives), and $F N$ (False Negatives) are the values from the confusion matrix.


#figure(
  image(".imgs/svm_confusion_matrix.svg"),
  caption: [Confusion matrix for the SVM model],
) <mlp-confusion-matrix>


= Results
== SP Motifs
From all SP-Endowed sequences, a motif logo was generated using the context around the cleavage site (13 bp upstream to 2bp downstream) @sp_motif.

#place(top + center, scope: "parent", float: true, [#figure(
    image(".imgs/positive_logo.svg"),

    caption: [Motif logo for the context (-13bp to 2bp) surrounding the cleavage site],
  )<sp_motif>,])

For the Von-Hejne model the computed optimal threshold was $9.25$, confusion matrix for the test set is reported in @vh-confusion-matrix.

#figure(
  image(".imgs/von_heijne_confusion_matrix.svg"),
  caption: [Confusion matrix for the Von-Heijne model],
) <vh-confusion-matrix>



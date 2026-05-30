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

= Materials and Methods
== Data Collection
Protein sequences were queried from UniprotKB, and collected into 2 distinct sets, positive and negative.
For Positive sequences, we selected sequences having experimental evidence for the presence of a signal peptide (`ft_signal_exp:*`)
For the negative set, proteins not annotated with a signal peptide and endowed with experimental evidence for subcellular localization in non-secretory comparments, namely nucleus (`cc_scl_term_exp:SL-0191`), peroxisome (`cc_scl_term_exp:SL-0204`), cell membrane / plasma membrane (`cc_scl_term_exp:SL-0039`), cytosol (`cc_scl_term_exp:SL-0091`), plastid (`cc_scl_term_exp:SL-0209`) or mitochondrion (`cc_scl_term_exp:SL-0173`).

Both sets were also limited to the taxa of eukaroyot and areviewed entries, and non-fragment nature.

Further selection on the positive entries was made to select proteins which had informaton about the cleavage site and a Signal Peptide with a length of at least 14bp.

// Initial filtering: 2960 entries to 2942
// After clustering 1093
For the negative set, additional information was added to the data representing the presence of a transmembrane helix, whose  hydrophobicity profile mimicks the one of a signal peptide.

// 20806 total, 18301 w/o TM, 2505 with TM
// after clustering: 9028

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


The threshold was computed using 5-fold cross validation.

== Clustering
Clustering was performed using mmseqs2 [cit], using parameters `=c 0.4 --min-seq-id 0.3 --cov-mode 0 --cluster-mode 1`, giving clusters with less than 30\% pairwise identity.


= Results
For the Von-Hejne model the computed optimal threshold was $9.25$.
// Image for the confusion matrix.


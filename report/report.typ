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


== Clustering
Clustering was performed using mmseqs2 [cit], using parameters `=c 0.4 --min-seq-id 0.3 --cov-mode 0 --cluster-mode 1`, giving clusters with less than 30\% pairwise identity.

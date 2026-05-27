#set document(
  title: [
    A Comparison of Con Heijne, MLP, and SVM Models for the classification of Signal Peptides.
  ],
)

#title()

= Materials and Methods
== Data Collection
Protein sequences were queried from UniprotKB, and collected into 2 distinct sets, positive and negative.
For Positive sequences, we selected sequences having experimental evidence for the presence of a signal peptide (`ft_signal_exp:*`)
For the negative set, proteins not annotated with a signal peptide and endowed with experimental evidence for subcellular localization in non-secretory comparments, namely nucleus (`cc_scl_term_exp:SL-0191`), peroxisome (`cc_scl_term_exp:SL-0204`), cell membrane / plasma membrane (`cc_scl_term_exp:SL-0039`), cytosol (`cc_scl_term_exp:SL-0091`), plastid (`cc_scl_term_exp:SL-0209`) or mitochondrion (`cc_scl_term_exp:SL-0173`).

Both sets were also limited to the taxa of eukaroyot and areviewed entries, and non-fragment nature.

Further selection on the positive entries was made to select proteins which had informaton about the cleavage site and a Signal Peptide with a length of at least 14bp.

For the negative set, additional information was added to the data representing the presence of a transmembrane helix, whose  hydrophobicity profile mimicks the one of a signal peptide.

== Clustering
Clustering was performed using mmseqs2 [cit],

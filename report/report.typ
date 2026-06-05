#set page(
  columns: 2,
)

#set text(
  size: 11pt,
)
#set heading(
  numbering: "1.",
)

#set par(
  leading: 0.56em,
  justify: true,
)

#place(
  top + left,
  scope: "parent",
  float: true,
  [
    #line(length: 100%)
    #v(0.5em)
    #text(1.6em, weight: "bold")[
      A comparative analysis of Von-Heijne, SVM and Perceptron Models for the detecton of Signal Peptides.

    ]
    #v(0.25em)
    #text(1.2em)[Fabrizio Guidotti#super("1")]

    // Affiliation
    #v(0.25em)
    #text(1em)[#super("1")Department of Pharmacy and Biotechnology, University of Bologna, Italy]

    #v(0.5em)
    #text(weight: "bold")[Abstract:] #linebreak()
    #text(weight: "bold")[Motivation:] Signal peptides enable proteins to enter the secretory pathway. Given their importance and the high cost of in-vivo assays, accurate in-silico SP detection models have taken on an important role in bioinformatics. #linebreak()
    #text(weight: "bold")[Results:] With this paper, we provide two models for the detection of Signal Peptides from sequence information. The first model, based on the Von-Heijne algorithm, produces a relatively accurate model, with an MCC of $0.66$. The second model is based on Support Vector Machines and achieves an higher MCC of $0.86$ and $f_1$ score of $0.87$, supporting an overall stronger model. Given the complexity of the full SVM model, we also trained a reduced model, which trades some performance (MCC of $0.73$, and $f_1$ of $0.76$), for a significantly lower number of features necessary, lowering the amount of preprocessing needed and the overall complexity of the final model.#linebreak()
    #text(weight: "bold")[Code Availability]: The full source code for this paper is available at #link("https://github.com/FabrizioG202/lab2").

    #v(0.5em)
    #line(length: 100%)
  ],
)

= Introduction
For a protein to enter the secretory pathway, in both eukaryotic and prokaryotic cells, it must be endowed with a specific target signal. Often, this signal takes the shape of a short sequences located at the N-terminus of proteins @von_Heijne_1990. The signal peptide (SP) has a distinct three-domain structure, depicted in @sp_structure, with a positively charged N-terminal region (n-region), a central hydrophobic region (H-region) and a more polar C-terminal region (C-region) containing the cleavage site @von_Heijne_1990. Given their importance in many aspects of cell biology, the accurate detection of SPs is a crucial task in bioinformatics, which has been tackled before with a variety of approaches, including machine learning (SVMs @Nugent_Jones_2009, @Petersen_Brunak_von_Heijne_Nielsen_2011 and Bayesian networks @Reynolds_Käll_Riffle_Bilmes_Noble_2008, Hidden Markov Models @Käll_Krogh_Sonnhammer_2004), and, more recently, deep learning @Savojardo_Martelli_Fariselli_Casadio_2018.

#figure(
  image(".imgs/sp_structure.png"),
  caption: [Signal peptide structure, adapted from @von_Heijne_1990],
) <sp_structure>


With this work, we propose a comparative analysis of two different approaches for the detection of SPs, namely the Von-Heijne model a Support Vector Machine (SVM) based method.


= Materials and Methods
== Data Collection
Reviewed non-fragment protein sequences from the eukaryota taxa were queried from UniprotKB @UniProt_Consortium_2018, and collected into 2 distinct sets, positive and negative.
For Positive sequences, we selected sequences having experimental evidence for the presence of a signal peptide (`ft_signal_exp:*`)
For the negative set, proteins not annotated with a signal peptide and endowed with experimental evidence for subcellular localization in non-secretory comparments, namely nucleus (`cc_scl_term_exp:SL-0191`), peroxisome (`cc_scl_term_exp:SL-0204`), cell membrane / plasma membrane (`cc_scl_term_exp:SL-0039`), cytosol (`cc_scl_term_exp:SL-0091`), plastid (`cc_scl_term_exp:SL-0209`) or mitochondrion (`cc_scl_term_exp:SL-0173`). For sequences in the negative set, we gathered information about the presence of a transmembrane (TM) helix, whose  hydrophobicity profile mimicks the one of a signal peptide.

Further selection on the positive entries was carried out in order to select proteins which had informaton about the cleavage site and a Signal Peptide with a length of at least 14bp.

To perform an initial screening for possible differences between the two sets, we plotted the distribution of sequence lengths in the positive and negative datasets (@length-distribution), which showed similar trends.

#figure(
  image(".imgs/length_distribution.svg"),
  caption: [Length distribution of positive and negative sequences],
) <length-distribution>

To ensure equal distribution among the kindoms, the composition of each set was considered (@kingdom-distribution).

#figure(
  image(".imgs/kingdom_distribution.svg"),
  caption: [Kingdom distribution of positive and negative datasets],
) <kingdom-distribution>

The distribution of cleavage-site positions, corresponding to the length of the cleaved N-terminal region in positive sequences, is shown in @cleaved-region-length-distribution. The distribution was used to quickly scan for outliers and to crudely confirm the biological plausibility of the information.

#figure(
  image(".imgs/cleaved_region_length_distribution.svg"),
  caption: [Length distribution of the cleaved region in positive sequences],
) <cleaved-region-length-distribution>

== Clustering

Sequences were clustered using mmseqs2 @Steinegger_Söding_2017, using parameters `--min-seq-id 0.3, cluster-mode 1 -cov-mode 0 -c 0.4`, in order to cluster sequences with a similarity of at least 30\% and an alignment spanning at least 40\% of the total sequence length. The number of sequences in the datasets is reported in @seq-counts.

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

The optimal threshold for the dataset was computed using 5-fold cross validation. For each fold, the data was split into validation and test sets. On the validation set, we computed precision and recall values across all possible threshold values, then selected the threshold that maximized the F1 score (see @evaluation-and-comparison). This optimal threshold was then evaluated on the corresponding test set to compute performance metrics. The final threshold reported is the average of the optimal thresholds across all 5 folds.

== Feature Extraction
To build a model that can deal with sequences, they must be encoded into numerical representation. There are many possibility of extracting meaningful numerical features from a sequence, for the SVM and MLP model, we employed 2 kinds of features:
- Aminoacid composition up to a cutoff $k$ (`n_term_composition`) and from $k$ to an hardcoded limit $N$ (set to $90$) (`c_term_composition`)
- Protein Scales computed on the first $N$ aminoacids with a window size of $5$ (@feature-example). An Hydrophobicity scale (@Kyte_Doolittle_1982), alpha-helix tendency scale (@Chou_Fasman_1978), transmembrane tendency (@Zhao_London_2009), bulkiness (@Zimmerman_Eliezer_Simha_1968), and polarity (@Zimmerman_Eliezer_Simha_1968) were used.

For this experiment, $k$ was chosen to be the average position of the cut site, roughly 22 aminoacids from the N-Terminus.

#figure(
  image(".imgs/feature_extraction_example.svg"),
  caption: [Example of feature extraction for a given sequence. Protein scales are computed on the first $N$ aminoacids with a window size of 5],
) <feature-example>

== Support Vector Machine <svm-methods>
Support Vector Machines (SVM) are a machine learning methodology that can be used for two-group classification problems. SVMs separate classes using an hyperplane, defined by a small set of points called support vectors. Albeit trivial for linearly separable classes, vector machines can make use of a kernel function, which allows the mapping of input points from the original space to an higher-dimension space (feature space). In this new space, a linear decision boundary can be established, which, when mapped back to the original space, will produce a complex non-linear decision boundary, enabling correct classification of complex classes @Cortes_Vapnik_1995. A commonly used kernel function is Radial Basis Function (RBF), which is defined as follows:
$
  K(x, x') = exp(-gamma ||x - x'||^2)
$

The $gamma$ parameter defines how much influence a single training example has. A low $gamma$ value means that the influence of a single training example is far-reaching, while a high $gamma$ value means that the influence is limited to close neighbors. For both linear and RBF kernels, the $C$ parameter controls the trade-off between maximizing the margin and minimizing the classification error, with low C values leading to a widwe margin and a simpler model, and high C values leading to a narrower margin, with few training errors but a more complex model.

We used a grid search to find optimal hyperparameters, testing both RBF and linear kernels. For the RBF kernel, we evaluated C values of 0.1, 1, 2, 4, and 8, combined with gamma values of "scale", 0.001, 0.01, 0.1, 1, and 2. For the linear kernel, we tested C values of 0.1, 1, 2, 4, and 8. Each model was tested on a and evaluated on an 80\% training split. The model yielding the highest f1 score was found to have a `rbf` kernel, with $C=2$ and $gamma=0.01$, with a mean f1 score across the 5 folds ran on the training set of $0.87$.

== Evaluation and comparison <evaluation-and-comparison>
Rather than relying on accuracy alone, since the dataset is greatly imbalanced, to evaluate and compare the performance of the models, we computed several metrics on the test set. The metrics were calculated based on the following formulas:

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


= Results
== SP Motifs
From all SP-Endowed sequences, a motif logo (@sp_motif) was generated using the context around the cleavage site (13 aa upstream to 2 aa downstream) .

#place(top + center, scope: "parent", float: true, [#figure(
    image(".imgs/positive_logo.svg"),

    caption: [Motif logo for the context (-13aa to 2aa) surrounding the cleavage site],
  )<sp_motif>,])

== Von Heijne Model

For the Von-Hejne model the computed optimal threshold was $9.25$, confusion matrix computed on the test set is reported in @vh-confusion-matrix. The model achieved an accuracy of $0.9378$, precision of $0.7474$, recall of $0.6532$, F1 score of $0.6971$, and MCC of $0.6645$. The model thus reaches an high accuracy, but performance on the positive class is limited, detecting only around 65% of actual positives.

#figure(
  image(".imgs/von_heijne_confusion_matrix.svg"),
  caption: [Confusion matrix for the Von-Heijne model],
) <vh-confusion-matrix>

== Support Vector Machine
The model, as chosen in the hyperparameter tuning phase, is a SVM with RBF kernel, $C=2$ and $gamma=0.01$. The confusion matrix computed on the test set is reported in @svm-confusion-matrix. The model achieved an accuracy of $0.9743$, precision of $0.9034$, recall of $0.8539$, F1 score of $0.8779$, and MCC of $0.8640$. The SVM model thus outperforms the Von-Heijne model across all metrics, showing particularly strong improvement in recall (85.39% vs 65.32%), indicating better detection of signal peptides. This means that the SVM model is both more reliable when predicting positives and much better at detecring actual positive cases. It is worth noting that this increase in performance mirrors a great increase in complexity of the underlying method, which is structurally more complex than a simple PSWM, and also requires far greater preprocessing to compute the features. To mitigate this, we tried producing a reduced model, using only the top 5 features selected using permutation importance.

#figure(
  image(".imgs/reduced_svm_confusion_matrix.svg"),
  caption: [Confusion matrix for the reduced SVM model],
) <reduced-svm-confusion-matrix>


=== Feature Importance <feature-importance>
In order to understand which features contribute the most to the ability of our model to classify sequences, we employed two methods. As a first test, we assesed permutation importance, where values of each feature are randomly shuffled in order to break their relationship with the target label, using the drop in the classifier accuracy as a proxy for the importance of the feature. As a second test, we trained a Random Forest Classifier with $500$ classifiers and exploited the model's own feature importance, which is computed as a result of the training process. We also evaluated the performance of the resulting random forest classifier, which itself proved to be a good model, with an accuracy of $0.9580$, precision of $0.9036$, recall of $0.6849$, F1 score of $0.7792$, and MCC of $0.7654$.

The top 5 features identified by each method are shown in @permutation-importance and @rf-importance. Both methods consistently identified `n_terminal_comp_L` (leucine composition in the N-terminal region) and `tm_max` (maximum transmembrane tendency) as important features, even if with different ranks.

#figure(
  table(
    columns: 2,
    align: (left, center),
    [*Feature*], [*Importance Mean*],
    [`n_terminal_comp_L`], [0.068479],
    [`n_terminal_comp_A`], [0.010622],
    [`tm_max`], [0.008453],
    [`kd_max_pos`], [0.008225],
    [`n_terminal_comp_D`], [0.008195],
  ),
  caption: [Top 5 features by permutation importance],
)<permutation-importance>

#figure(
  table(
    columns: 2,
    align: (left, center),
    [*Feature*], [*Importance*],
    [`kd_max_pos`], [0.11984],
    [`tm_max_pos`], [0.105942],
    [`tm_max`], [0.091583],
    [`kd_max`], [0.078178],
    [`n_terminal_comp_L`], [0.06908],
  ),
  caption: [Top 5 features by random forest importance],
)<rf-importance>

=== Reduced SVM Model
We then trained a reduced model, using the combination of 5 most important features found by each method. The chosen features were:
- `kd_max`: Maximum hydrophobicity value
- `kd_max_pos`: Position of maximum hydrophobicity
- `n_terminal_comp_A`: Alanine composition in the N-terminal region
- `n_terminal_comp_D`: Aspartic acid composition in the N-terminal region
- `n_terminal_comp_L`: Leucine composition in the N-terminal region
- `tm_max`: Maximum transmembrane tendency
- `tm_max_pos`: Position of maximum transmembrane tendency

A model was selected using the same method used to select the full model, including the same pool of possible classificators. The reduced model achieved an accuracy of $0.9521$, precision of $0.8177$, recall of $0.7169$, F1 score of $0.7640$, and MCC of $0.7394$. The confusion matrix computed on the test set is reported in @reduced-svm-confusion-matrix. While this represents a decrease in performance compared to the full SVM model (F1 score of 0.8779 vs 0.7640), the reduced model still outperforms the Von-Heijne model (F1 score of 0.6971) while using only 7 features instead of the full feature set.


#figure(
  image(".imgs/svm_confusion_matrix.svg"),
  caption: [Confusion matrix for the SVM model],
) <svm-confusion-matrix>


= Code Availability
The code used for this analysis is available at: #link("https://github.com/FabrizioG202/lab2")



#bibliography("bib.bib")

import logomaker
import matplotlib.pyplot as plt

# ( These are just randomly generated color, since chemistry based colors do not seem to work)
_AA_COLORS = {
    "A": "#8DD3C7",
    "R": "#FB8072",
    "N": "#80B1D3",
    "D": "#FDB462",
    "C": "#B3DE69",
    "Q": "#FCCDE5",
    "E": "#BC80BD",
    "G": "#CCEBC5",
    "H": "#FFED6F",
    "I": "#1B9E77",
    "L": "#D95F02",
    "K": "#7570B3",
    "M": "#E7298A",
    "F": "#66A61E",
    "P": "#E6AB02",
    "S": "#A6761D",
    "T": "#666666",
    "W": "#A6CEE3",
    "Y": "#B2DF8A",
    "V": "#FB9A99",
}


def generate_logo(
    sequences: list[str], left_bound: int, right_bound: int
) -> plt.Figure:
    # Compute the counts for the logo.
    counts_matrix = logomaker.alignment_to_matrix(
        sequences,
        to_type="information",
        characters_to_ignore="X",
    )

    fig, ax = plt.subplots(figsize=(12, 4))

    # Set x axis as relative to the cleavage site.
    ax.set_xticks(range(left_bound + right_bound))
    ax.set_xticklabels(
        [str(i) for i in range(-left_bound, 0)]
        + [str(i) for i in range(0, right_bound)],
    )

    logomaker.Logo(
        counts_matrix,
        shade_below=0.5,
        fade_below=0.5,
        font_name="Arial Rounded MT Bold",
        ax=ax,
        color_scheme=_AA_COLORS,
    )

    return fig

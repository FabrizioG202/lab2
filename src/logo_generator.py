import logomaker
import matplotlib.pyplot as plt


def generate_logo(
    sequences: list[str], left_bound: int, right_bound: int
) -> plt.Figure:
    # Display the motif as a sequence logo
    counts_matrix = logomaker.alignment_to_matrix(
        sequences,
        to_type="information",
        characters_to_ignore="X",
    )

    fig, ax = plt.subplots(figsize=(12, 4))
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
    )

    return fig

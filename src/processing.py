import polars as pl


def split_df(
    df: pl.DataFrame, fraction: float, seed: int
) -> tuple[pl.DataFrame, pl.DataFrame]:
    # create a shuffled version of the dataframe
    shuffled = df.sample(fraction=1.0, shuffle=True, seed=seed)

    # split it into 2
    count = int(len(df) * fraction)
    return shuffled[:count], shuffled[count:]

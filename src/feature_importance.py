import plotly.graph_objects as go
import polars as pl
import sklearn.base


def draw_feature_importance(
    svm_model: sklearn.base.BaseEstimator,
    x_train: pl.DataFrame,
    y_train: pl.Series,
    x_dataframe: pl.DataFrame,
) -> go.Figure:

    # First, compute importance using permurations
    result = sklearn.inspection.permutation_importance(
        svm_model,
        x_train,
        y_train,
        scoring="f1_macro",
        n_repeats=10,
        random_state=42,
        n_jobs=-1,
    )

    permutation_importance = pl.DataFrame(
        {
            "feature": x_dataframe.columns,
            "importance_mean": result.importances_mean,
            "importance_std": result.importances_std,
        }
    ).sort("importance_mean", descending=True)

    # Now, Compute the importance using random forest feature importance
    random_forest = sklearn.ensemble.RandomForestClassifier(
        n_estimators=500, random_state=42, n_jobs=-1, class_weight="balanced"
    )

    # Fit the model
    _ = random_forest.fit(x_train, y_train)

    random_forest_importance = pl.DataFrame(
        {
            "feature": X_df.columns,
            "importance": random_forest.feature_importances_,
        }
    ).sort("importance", descending=True)

    # plot the pi and rf importnace as side by side horizontal bar plots using plotly
    fig = go.Figure(
        data=[
            go.Bar(
                x=permutation_importance["importance_mean"],
                y=permutation_importance["feature"],
                orientation="h",
                name="Permutation Importance",
                error_x=dict(
                    type="data", array=permutation_importance["importance_std"]
                ),
            ),
            go.Bar(
                x=random_forest_importance["importance"],
                y=random_forest_importance["feature"],
                orientation="h",
                name="Random Forest Importance",
            ),
        ]
    )

    fig.update_layout(
        title="Feature Importance",
        xaxis_title="Importance",
        yaxis_title="Feature",
        yaxis={
            "categoryorder": "total ascending",
            "tickfont": {"size": 10},
            "tickmode": "linear",  # So that all the labels are shown, even if they overlap
        },
        barmode="group",
    )

    return fig

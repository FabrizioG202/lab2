from dataclasses import dataclass

import plotly.graph_objects as go
import polars as pl
import sklearn.base
import sklearn.ensemble
import sklearn.inspection


@dataclass
class FeatureImportance:
    feature_df: pl.DataFrame
    x_train: object
    y_train: object
    best_model: sklearn.base.BaseEstimator
    label_column: str = "label"
    id_columns: tuple[str, ...] = ("accession",)

    @property
    def feature_columns(self) -> list[str]:
        return [
            col
            for col in self.feature_df.columns
            if col not in {self.label_column, *self.id_columns}
        ]

    def get_feature_matrix(self) -> pl.DataFrame:
        return self.feature_df.select(self.feature_columns)

    def get_labels(self):
        return self.feature_df[self.label_column].to_numpy()

    def compute_permutation_importance(
        self,
        *,
        scoring: str = "f1_macro",
        n_repeats: int = 1,
        random_state: int = 42,
    ) -> pl.DataFrame:
        result = sklearn.inspection.permutation_importance(
            self.best_model,
            self.x_train,
            self.y_train,
            scoring=scoring,
            n_repeats=n_repeats,
            random_state=random_state,
            n_jobs=-1,
        )

        return pl.DataFrame(
            {
                "feature": self.feature_columns,
                "importance_mean": result.importances_mean,
                "importance_std": result.importances_std,
            }
        ).sort("importance_mean", descending=True)

    def compute_random_forest_importance(
        self,
        *,
        n_estimators: int = 500,
        random_state: int = 42,
    ) -> tuple[pl.DataFrame, sklearn.base.BaseEstimator]:
        random_forest = sklearn.ensemble.RandomForestClassifier(
            n_estimators=n_estimators,
            random_state=random_state,
            n_jobs=-1,
            class_weight="balanced",
        )
        _ = random_forest.fit(self.x_train, self.y_train)

        return (
            pl.DataFrame(
                {
                    "feature": self.feature_columns,
                    "importance": random_forest.feature_importances_,
                }
            ).sort("importance", descending=True),
            random_forest,
        )

    def draw_feature_importance(
        self,
        permutation_importance: pl.DataFrame,
        random_forest_importance: pl.DataFrame,
    ) -> go.Figure:
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
            title="Feature Importance (Permutation Importance and Random Forest)",
            xaxis_title="Importance",
            yaxis_title="Feature",
            yaxis={
                "categoryorder": "total ascending",
                "tickfont": {"size": 10},
                "tickmode": "linear",
            },
            barmode="group",
        )

        return fig

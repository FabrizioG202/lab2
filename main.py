import re

from src.data_collection import DataCollector
from src.utils import AdditionalProtParamData

collector = DataCollector(
    positive_query="""
    (
        (taxonomy_id:2759)
        AND (reviewed:true)
        AND (ft_signal_exp:*)
        AND (fragment:false)
        AND (length:[40 TO *])
        AND (existence:1)
    )
    """,
    negative_query="""
    (
        (reviewed:true)
        AND (fragment:false)
        AND (taxonomy_id:2759)
        AND (length:[40 TO *])
        AND (existence:1)
        AND NOT (ft_signal:*)
        AND (
            (cc_scl_term_exp:SL-0191)
            OR (cc_scl_term_exp:SL-0204)
            OR (cc_scl_term_exp:SL-0039)
            OR (cc_scl_term_exp:SL-0091)
            OR (cc_scl_term_exp:SL-0209)
            OR (cc_scl_term_exp:SL-0173)
        )
    )
    """,
)


# Setup the working director for the project
# (creates .data, and .imgs)
collector.setup_wd()


# Get the positive examples and cluster them
positive = collector.cluster_df(collector.get_positive_examples())

# Get the negative examples and cluster them
negative = collector.cluster_df(collector.get_negative_examples())

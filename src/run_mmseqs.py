import subprocess


def run_mmseqs(
    input_path: str,
    output_path: str,
    /,
    tmp_path: str = "tmp",
    min_seq_id: float = 0.3,
    cluster_mode: int = 1,
    cov_mode: int = 0,
    c: float = 0.4,
) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            "mmseqs",
            "easy-cluster",
            input_path,
            output_path,
            tmp_path,
            "--min-seq-id",
            str(min_seq_id),
            "-c",
            str(c),
            "--cov-mode",
            str(cov_mode),
            "--cluster-mode",
            str(cluster_mode),
            # "--remove-tmp-files",
            # "true",
        ],
        check=True,
        capture_output=True,
        shell=True,
    )

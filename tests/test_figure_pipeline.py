import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "roary_pangenome_figures.py"
INPUT = REPO / "data" / "core_genes_reps.csv"


class FigurePipelineTest(unittest.TestCase):
    def test_representative_panel_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "figures"
            env = os.environ.copy()
            env["MPLCONFIGDIR"] = str(Path(temp) / "matplotlib")
            subprocess.run(
                [sys.executable, str(SCRIPT), "--input", str(INPUT), "--outdir", str(output), "--seed", "42"],
                check=True, cwd=REPO, env=env,
            )
            for name in [
                "01_prevalence_distribution.png", "02_presence_absence_heatmap.png",
                "03_genome_gene_content.png", "04_prevalence_copy_number.png",
                "05_gene_content_pcoa.png", "gene_family_prevalence.tsv",
                "genome_content.tsv", "run_manifest.json", "report.html",
            ]:
                self.assertTrue((output / name).is_file(), name)
            self.assertFalse((output / "06_permutation_accumulation.png").exists())


if __name__ == "__main__":
    unittest.main()

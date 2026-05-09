from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]

REQUIRED = [
    "README.md",
    "baseline/mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_20260508.xlsx",
    "baseline/mymodel_CGA009_publishable_v1.13_lps_acp_acyltransferase_for_memote.xml",
    "baseline/CLAUDE_HANDOFF_V1.13_20260509.md",
    "baseline/v1.13_publication_summary_metrics.tsv",
    "baseline/cobra_matlab_qc_summary.tsv",
    "baseline/cga009_phenotype_profile_results.tsv",
    "evidence/literature_evidence.tsv",
    "evidence/LITERATURE_READING_PROTOCOL.md",
    "evidence/phenotype_data_requirements.tsv",
    "next_agent/NEXT_AGENT_PROMPT.md",
    "skills/cga009-top-journal-modeling/SKILL.md",
]


def main() -> int:
    missing = [path for path in REQUIRED if not (ROOT / path).exists()]
    if missing:
        print("Missing required files:")
        for path in missing:
            print(f"- {path}")
        return 1
    print(f"OK: {len(REQUIRED)} required files present under {ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

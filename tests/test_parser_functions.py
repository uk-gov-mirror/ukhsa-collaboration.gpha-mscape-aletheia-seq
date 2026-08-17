import datetime

import pandas as pd
import pytest

from AletheiaSeq import aletheiaseq_parse

"""
FIXTURES
"""


@pytest.fixture
def valid_yaml():
    yml = {
        "SchemaDetails": {"name": "Streptococcus", "type": "speciation", "version": "v0.1"},
        "BlastDetails": {
            "seqs": ["ply", "lytA", "SP2020", "gmuR", "wzg", "kdpE", "scpC", "cfb", "fbpS"],
            "db_name": "Streptococcus_speciation",
            "db_type": "nucleotide",
            "last_update": datetime.date(2026, 8, 6),
            "outfmt": [
                "qseqid",
                "sseqid",
                "qlen",
                "slen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "length",
                "pident",
            ],
            "filters": ["pident >= 80", "evalue <= 1e-10", "length >= (slen * 0.8)"],
            "min_hits": 3,
        },
        "SchemaGroups": [
            {
                "name": "test group",
                "defining_loci": ["locus1", "locus2", "locus3"],
                "required_loci": ["locus1"],
                "minimum_loci_required": 2,
            },
        ],
    }
    return yml


@pytest.fixture
def parser(valid_yaml):
    p = aletheiaseq_parse.AletheiaSeqParser.__new__(aletheiaseq_parse.AletheiaSeqParser)

    p.yml = valid_yaml

    p.loci_results = {}
    return p


@pytest.fixture
def blast_df():
    blast_df = pd.DataFrame(
        {
            "qseqid": ["pass_all", "fail_length", "fail_pident", "fail_evalue"],
            "sseqid": ["locus1", "locus1", "locus1", "locus1"],
            "slen": [160, 160, 160, 160],
            "pident": [95, 80, 65, 85],
            "evalue": [0.0, 0.0, 0.0, 3],
            "length": [150, 10, 159, 160],
        }
    )
    return blast_df


@pytest.fixture
def blast_summary_pass():
    blast_summary = pd.DataFrame(
        {
            "sseqid": ["locus1", "locus3"],
            "total_hits": [2, 2],
            "passed_hits": [2, 1],
        }
    )
    return blast_summary


"""
Test parser
"""


def test_filter_blast_applies_filters(parser, blast_df):
    parser.blast_df = blast_df

    parser._AletheiaSeqParser__filter_blast()

    assert parser.blast_df["passed_filters"].tolist() == [True, False, False, False]


def test_filter_blast_no_filters_marks_all_pass(parser, blast_df):
    parser.yml["BlastDetails"]["filters"] = []
    parser.blast_df = blast_df

    parser._AletheiaSeqParser__filter_blast()

    assert parser.blast_df["passed_filters"].all()


def test_report_negative_populates_all_groups(parser):
    parser._AletheiaSeqParser__report_negative()

    result = parser.loci_results["test group"]

    assert result["all_req"] is False
    assert result["min_loci"] == 0
    assert result["min_loci_pass"] is False


def test_report_negative_loci(parser):
    parser._AletheiaSeqParser__report_negative()

    df = parser.loci_results["test group"]["results_df"]

    assert set(df["sseqid"]) == {"locus1", "locus2", "locus3"}
    assert df.loc[(df.sseqid == "locus1"), "required"].item() is True
    assert df.loc[(df.sseqid == "locus2"), "required"].item() is False
    assert df.loc[(df.sseqid == "locus3"), "required"].item() is False


def test_summarise_groups_adds_missing_defining_loci(parser, blast_summary_pass):
    parser.blast_summary = blast_summary_pass
    parser._AletheiaSeqParser__summarise_groups()

    df = parser.loci_results["test group"]["results_df"]

    assert set(df["sseqid"]) == {
        "locus1",
        "locus2",
        "locus3",
    }

    locus2 = df[df.sseqid == "locus2"].iloc[0]

    assert locus2["total_hits"] == 0
    assert locus2["passed_hits"] == 0


def test_summarise_groups_calculates_threshold(parser, blast_summary_pass):
    parser.blast_summary = blast_summary_pass

    parser._AletheiaSeqParser__summarise_groups()

    df = parser.loci_results["test group"]["results_df"]

    assert df.loc[(df.sseqid == "locus1"), "above_min_threshold"].item() is True
    assert df.loc[(df.sseqid == "locus3"), "above_min_threshold"].item() is False


def test_summarise_groups_all_required_loci_present(parser, blast_summary_pass):
    parser.blast_summary = blast_summary_pass

    parser._AletheiaSeqParser__summarise_groups()

    result = parser.loci_results["test group"]

    assert result["all_req"] is True


def test_summarise_groups_required_locus_missing(parser):
    parser.blast_summary = pd.DataFrame(
        {
            "sseqid": ["locus1", "locus3"],
            "total_hits": [2, 2],
            "passed_hits": [1, 2],
        }
    )

    parser._AletheiaSeqParser__summarise_groups()

    result = parser.loci_results["test group"]

    assert result["all_req"] is False


def test_summarise_groups_minimum_loci_passes(parser):
    parser.blast_summary = pd.DataFrame(
        {
            "sseqid": ["locus1", "locus2"],
            "total_hits": [2, 2],
            "passed_hits": [2, 2],
        }
    )

    parser._AletheiaSeqParser__summarise_groups()

    result = parser.loci_results["test group"]

    assert result["min_loci"] == 2
    assert result["min_loci_pass"] is True


def test_summarise_groups_minimum_loci_fails(parser):
    parser.blast_summary = pd.DataFrame(
        {
            "sseqid": ["locus1"],
            "total_hits": [2],
            "passed_hits": [2],
        }
    )

    parser._AletheiaSeqParser__summarise_groups()

    result = parser.loci_results["test group"]

    assert result["min_loci"] == 1
    assert result["min_loci_pass"] is False


def test_empty_required_loci_always_passes(parser, monkeypatch, blast_df):
    parser.yml["SchemaGroups"][0]["required_loci"] = []

    monkeypatch.setattr(
        parser,
        "_AletheiaSeqParser__load_blast",
        lambda _: setattr(parser, "blast_df", blast_df),
    )

    parser.process_blast("dummy.tsv")

    assert parser.loci_results["test group"]["all_req"] is True


def test_process_blast_empty_file_generates_negative_results(parser, monkeypatch):
    parser.blast_df = pd.DataFrame()

    monkeypatch.setattr(
        parser,
        "_AletheiaSeqParser__load_blast",
        lambda _: None,
    )

    parser.process_blast("dummy.tsv")

    result = parser.loci_results["test group"]

    assert result["all_req"] is False
    assert result["min_loci"] == 0
    assert result["min_loci_pass"] is False


def test_process_blast_populates_loci_results(parser, blast_df, monkeypatch):
    parser.blast_df = blast_df

    monkeypatch.setattr(
        parser,
        "_AletheiaSeqParser__load_blast",
        lambda _: None,
    )

    parser.process_blast("dummy.tsv")

    result = parser.loci_results["test group"]
    df = result["results_df"]

    assert result["all_req"] is False
    assert result["min_loci"] == 0
    assert result["min_loci_pass"] is False
    assert set(df["sseqid"]) == {
        "locus1",
        "locus2",
        "locus3",
    }


def test_process_blast_full_workflow(parser, monkeypatch, blast_df):
    monkeypatch.setattr(
        parser,
        "_AletheiaSeqParser__load_blast",
        lambda _: setattr(parser, "blast_df", blast_df),
    )

    parser.process_blast("dummy.tsv")

    result = parser.loci_results["test group"]

    assert result["all_req"] is False
    assert result["min_loci"] == 0

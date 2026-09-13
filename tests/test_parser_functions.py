import datetime

import pandas as pd
import pytest

from AletheiaSeq import aletheiaseq_parse, aletheiaseq_validate_yml

"""
FIXTURES
"""


@pytest.fixture
def valid_yaml():
    yml = {
        "SchemaDetails": {"name": "Test", "schema_type": "speciation", "version": "0.1.0"},
        "BlastDetails": {
            "seqs": [
                "locus1",
                "locus2",
                "locus3",
            ],
            "db_name": "Test_speciation",
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

    p.yml = aletheiaseq_validate_yml.SchemaConfig.model_validate(valid_yaml)

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
            "total_hits": [5, 3],
            "passed_hits": [3, 2],
        }
    )
    return blast_summary


@pytest.fixture
def schema_group_summary():
    locus1 = aletheiaseq_parse.LocusResult(
        locus_name="locus1", total_hits=5, passed_hits=3, above_min_threshold=True, required=True
    )
    locus2 = aletheiaseq_parse.LocusResult(
        locus_name="locus2", total_hits=5, passed_hits=3, above_min_threshold=True, required=False
    )
    locus3 = aletheiaseq_parse.LocusResult(
        locus_name="locus3", total_hits=5, passed_hits=1, above_min_threshold=False, required=False
    )
    summary = aletheiaseq_parse.SchemaGroupSummary(
        loci=[locus1, locus2, locus3], all_req=True, detected_loci_count=2, min_loci_pass=True
    )

    return summary


"""
Test parser
"""


def test_blast_load_column_mismatch(blast_df, tmp_path, parser):
    with pytest.raises(ValueError, match="but YAML outfmt specifies"):
        tmp_fp = tmp_path / "test_blast_out.tsv"
        blast_df.to_csv(tmp_fp, sep="\t", index=False)
        parser._AletheiaSeqParser__load_blast(tmp_fp)


def test_blast_load_successful(tmp_path, parser):
    tmp_fp = tmp_path / "test_blast_out.tsv"
    tmp_fp.write_text(
        "qseqid1\tsseqid1\t100\t150\t10\t100\t50\t140\t0\t0\t90\t85\nqseqid2\tsseqid2\t100\t150\t10\t100\t50\t140\t0\t0\t90\t85\n"
    )
    parser._AletheiaSeqParser__load_blast(tmp_fp)

    assert parser.blast_df.shape == (2, 12)


def test_blast_load_successful_empty_df(tmp_path, parser):
    tmp_fp = tmp_path / "test_blast_out.tsv"
    tmp_fp.write_text("\n")
    parser._AletheiaSeqParser__load_blast(tmp_fp)

    assert parser.blast_df.shape == (0, 12)


def test_filter_blast_applies_filters(parser, blast_df):
    parser.blast_df = blast_df

    parser._AletheiaSeqParser__filter_blast()

    assert parser.blast_df["passed_filters"].tolist() == [True, False, False, False]


def test_filter_blast_no_filters_marks_all_pass(parser, blast_df):
    parser.yml.BlastDetails.filters = []
    parser.blast_df = blast_df

    parser._AletheiaSeqParser__filter_blast()

    assert parser.blast_df["passed_filters"].all()


def test_report_negative_populates_all_groups(parser):
    parser.blast_df = pd.DataFrame(columns=parser.yml.BlastDetails.outfmt)
    parser._AletheiaSeqParser__report_negative()

    result = parser.loci_results["test group"]

    assert result.all_req is False
    assert result.detected_loci_count == 0
    assert result.min_loci_pass is False


def test_report_negative_loci(parser):
    parser.blast_df = pd.DataFrame(columns=parser.yml.BlastDetails.outfmt)
    parser._AletheiaSeqParser__report_negative()

    parser._AletheiaSeqParser__report_negative()

    df = pd.DataFrame(parser.loci_results["test group"].loci)

    assert set(df["locus_name"]) == {"locus1", "locus2", "locus3"}
    assert df.loc[(df.locus_name == "locus1"), "required"].item() is True
    assert df.loc[(df.locus_name == "locus2"), "required"].item() is False
    assert df.loc[(df.locus_name == "locus3"), "required"].item() is False


def test_summarise_groups_adds_missing_defining_loci(parser, blast_summary_pass):
    blast_summary = blast_summary_pass
    parser._AletheiaSeqParser__summarise_groups(blast_summary)

    df = pd.DataFrame(parser.loci_results["test group"].loci)

    assert set(df["locus_name"]) == {
        "locus1",
        "locus2",
        "locus3",
    }

    locus2 = df[df.locus_name == "locus2"].iloc[0]

    assert locus2["total_hits"] == 0
    assert locus2["passed_hits"] == 0


def test_summarise_groups_calculates_threshold(parser, blast_summary_pass):
    parser._AletheiaSeqParser__summarise_groups(blast_summary_pass)

    df = pd.DataFrame(parser.loci_results["test group"].loci)

    assert df.loc[(df.locus_name == "locus1"), "above_min_threshold"].item() is True
    assert df.loc[(df.locus_name == "locus3"), "above_min_threshold"].item() is False


def test_summarise_groups_all_required_loci_present(parser, blast_summary_pass):
    parser._AletheiaSeqParser__summarise_groups(blast_summary_pass)

    result = parser.loci_results["test group"]

    assert result.all_req is True


def test_summarise_groups_required_locus_missing(parser):
    blast_summary = pd.DataFrame(
        {
            "sseqid": ["locus1", "locus3"],
            "total_hits": [2, 2],
            "passed_hits": [1, 2],
        }
    )

    parser._AletheiaSeqParser__summarise_groups(blast_summary)

    result = parser.loci_results["test group"]

    assert result.all_req is False


def test_summarise_groups_minimum_loci_passes(parser):
    blast_summary = pd.DataFrame(
        {
            "sseqid": ["locus1", "locus2"],
            "total_hits": [3, 3],
            "passed_hits": [3, 3],
        }
    )

    parser._AletheiaSeqParser__summarise_groups(blast_summary)

    result = parser.loci_results["test group"]

    assert result.detected_loci_count == 2
    assert result.min_loci_pass is True


def test_summarise_groups_minimum_loci_fails(parser):
    blast_summary = pd.DataFrame(
        {
            "sseqid": ["locus1"],
            "total_hits": [3],
            "passed_hits": [3],
        }
    )

    parser._AletheiaSeqParser__summarise_groups(blast_summary)

    result = parser.loci_results["test group"]

    assert result.detected_loci_count == 1
    assert result.min_loci_pass is False


def test_empty_required_loci_always_passes(parser, monkeypatch, blast_df):
    parser.yml.SchemaGroups[0].required_loci = []

    monkeypatch.setattr(
        parser,
        "_AletheiaSeqParser__load_blast",
        lambda _: setattr(parser, "blast_df", blast_df),
    )

    parser.process_blast("dummy.tsv")

    assert parser.loci_results["test group"].all_req is True


def test_process_blast_empty_file_generates_negative_results(parser, monkeypatch):
    parser.blast_df = pd.DataFrame()

    monkeypatch.setattr(
        parser,
        "_AletheiaSeqParser__load_blast",
        lambda _: None,
    )

    parser.process_blast("dummy.tsv")

    result = parser.loci_results["test group"]

    assert result.all_req is False
    assert result.detected_loci_count == 0
    assert result.min_loci_pass is False


def test_process_blast_populates_loci_results(parser, blast_df, monkeypatch):
    parser.blast_df = blast_df

    monkeypatch.setattr(
        parser,
        "_AletheiaSeqParser__load_blast",
        lambda _: None,
    )

    parser.process_blast("dummy.tsv")

    result = parser.loci_results["test group"]
    locus_names = [locus.locus_name for locus in parser.loci_results["test group"].loci]

    assert result.all_req is False
    assert result.detected_loci_count == 0
    assert result.min_loci_pass is False
    assert set(locus_names) == {
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

    assert result.all_req is False
    assert result.detected_loci_count == 0


def test_schema_group_summary_headline(schema_group_summary):
    assert schema_group_summary.headline("test") == "test:required-pass,minimum-pass"

    schema_group_summary.min_loci_pass = False
    assert schema_group_summary.headline("test") == "test:required-pass,minimum-fail"

    schema_group_summary.all_req = False
    assert schema_group_summary.headline("test") == "test:required-fail,minimum-fail"


def test_schema_group_summary_summary(schema_group_summary):
    assert (
        schema_group_summary.summary()
        == "Of the 3 defining loci, 2 were detected and passed relevant thresholds. All required loci were detected."
    )

    passed_required_locus = aletheiaseq_parse.LocusResult(
        locus_name="passed_required_locus", total_hits=5, passed_hits=3, above_min_threshold=True, required=True
    )
    failed_required_locus = aletheiaseq_parse.LocusResult(
        locus_name="passed_required_locus", total_hits=5, passed_hits=2, above_min_threshold=False, required=True
    )
    passed_nonrequired_locus = aletheiaseq_parse.LocusResult(
        locus_name="passed_required_locus", total_hits=5, passed_hits=3, above_min_threshold=True, required=False
    )

    schema_group_summary.loci = [passed_nonrequired_locus]
    assert (
        schema_group_summary.summary()
        == "Of the 1 defining loci, 2 were detected and passed relevant thresholds. No loci are defined as required."
    )

    schema_group_summary.loci = [passed_required_locus, passed_nonrequired_locus]
    schema_group_summary.detected_loci_count = 2
    assert (
        schema_group_summary.summary()
        == "Of the 2 defining loci, 2 were detected and passed relevant thresholds. All required loci were detected."
    )

    schema_group_summary.loci = [passed_required_locus, failed_required_locus, passed_nonrequired_locus]
    schema_group_summary.all_req = False
    assert (
        schema_group_summary.summary()
        == "Of the 3 defining loci, 2 were detected and passed relevant thresholds. Not all required loci were detected."
    )


def test_html_outputs(parser, schema_group_summary, tmp_path):
    aletheiaseq_parse.create_html_output(
        total_hits=15,
        passed_hits=7,
        yaml_details=parser.yml.SchemaDetails,
        methods_dict=parser.yml.summarise_methods(),
        groups_matched=["test"],
        loci_results={"test": schema_group_summary},
        sample_id="ID-1234",
        output_directory=tmp_path.as_posix(),
    )

    with pytest.warns(UserWarning, match="HTML output directory already exists"):
        aletheiaseq_parse.create_html_output(
            total_hits=15,
            passed_hits=7,
            yaml_details=parser.yml.SchemaDetails,
            methods_dict=parser.yml.summarise_methods(),
            groups_matched=["test"],
            loci_results={"test": schema_group_summary},
            sample_id="ID-1234",
            output_directory=tmp_path.as_posix(),
        )

        assert (tmp_path / "html_summary_files/blast_headline.txt").exists()
        assert (tmp_path / "html_summary_files/test_locus_results.csv").exists()


def test_onyx_output(parser, schema_group_summary):
    onyx_analysis, exitcode = aletheiaseq_parse.create_analysis_fields(
        sample_id="ID-1234",
        yaml_details=parser.yml.SchemaDetails,
        yaml_methods_dict=parser.yml.summarise_methods(),
        headline_result=schema_group_summary.headline("test"),
        loci_results={"test": schema_group_summary},
        server="mscape",
    )

    assert exitcode == 0

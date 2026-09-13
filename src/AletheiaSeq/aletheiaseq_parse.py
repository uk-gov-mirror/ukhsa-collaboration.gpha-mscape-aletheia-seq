import logging
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

import pandas as pd
from onyx_analysis_helper import onyx_analysis_helper_functions as oa

from AletheiaSeq import aletheiaseq_validate_yml


@dataclass
class LocusResult:
    locus_name: str
    total_hits: int
    passed_hits: int
    above_min_threshold: bool
    required: bool


@dataclass
class SchemaGroupSummary:
    loci: list[LocusResult]
    all_req: bool
    detected_loci_count: int
    min_loci_pass: bool

    def headline(self, name):
        required = "pass" if self.all_req else "fail"

        min_loci = "pass" if self.min_loci_pass else "fail"

        return f"{name}:required-{required},minimum-{min_loci}"

    def summary(self):
        if self.all_req and any(locus.required for locus in self.loci):
            req = "All required loci were detected."
        elif self.all_req:
            req = "No loci are defined as required."
        else:
            req = "Not all required loci were detected."

        return f"Of the {len(self.loci)} defining loci, {self.detected_loci_count} were detected and passed relevant thresholds. {req}"


class AletheiaSeqParser:
    def __init__(self, yaml_file: str):
        self.yml: aletheiaseq_validate_yml.SchemaConfig = aletheiaseq_validate_yml.validate_yaml(yaml_file)
        self.loci_results: dict[str, SchemaGroupSummary] = {}

    def __load_blast(self, blast_out: str):
        columns = self.yml.BlastDetails.outfmt

        try:
            raw_df = pd.read_table(blast_out, sep="\t", header=None)

            expected = len(columns)
            actual = raw_df.shape[1]

            if actual != expected:
                raise ValueError(
                    f"BLAST output contains {actual} columns but YAML outfmt specifies {expected} columns."
                )
            else:
                raw_df.columns = columns
                self.blast_df = raw_df

                logging.info(
                    f"Loaded blast output file with column headers: {','.join(columns)}. Number of rows loaded: {self.blast_df.shape[0]}."
                )
        except pd.errors.EmptyDataError:
            self.blast_df = pd.DataFrame(columns=columns)

    def __filter_blast(self):
        filters = self.yml.BlastDetails.filters

        if len(filters) > 0:
            self.blast_df["passed_filters"] = self.blast_df.eval(" & ".join(filters))
        else:
            logging.info("No filters are present in the YAML file so all rows are considered QC pass")
            self.blast_df["passed_filters"] = True

    def __report_negative(self):
        self.blast_df["passed_filters"] = None
        for group in self.yml.SchemaGroups:
            results = []
            for loci in group.defining_loci:
                req = loci in group.required_loci
                locus = LocusResult(loci, 0, 0, False, req)
                results.append(locus)

            s = SchemaGroupSummary(loci=results, all_req=False, detected_loci_count=0, min_loci_pass=False)

            self.loci_results[group.name] = s

    def __summarise_groups(self, blast_summary: pd.DataFrame):
        summary_lookup = {row.sseqid: row for row in blast_summary.itertuples()}

        for group in self.yml.SchemaGroups:
            results = []
            required_hits = set()
            threshold_hits = set()

            for locus in group.defining_loci:
                row = summary_lookup.get(locus)

                total_hits = 0
                passed_hits = 0

                if row is not None:
                    total_hits = int(row.total_hits)
                    passed_hits = int(row.passed_hits)

                above_threshold = passed_hits >= self.yml.BlastDetails.min_hits

                required = locus in group.required_loci
                if above_threshold:
                    threshold_hits.add(locus)

                    if required:
                        required_hits.add(locus)

                results.append(
                    LocusResult(
                        locus_name=locus,
                        total_hits=total_hits,
                        passed_hits=passed_hits,
                        above_min_threshold=above_threshold,
                        required=required,
                    )
                )

            detected_loci_count = len(threshold_hits)

            summary = SchemaGroupSummary(
                loci=results,
                all_req=set(group.required_loci).issubset(required_hits),
                detected_loci_count=detected_loci_count,
                min_loci_pass=detected_loci_count >= group.minimum_loci_required,
            )

            self.loci_results[group.name] = summary

    def process_blast(self, blast_out: str):
        self.__load_blast(blast_out)

        if self.blast_df.shape[0] == 0:
            logging.info("No blast hits in output. Will report negative for all groups")
            self.__report_negative()
        else:
            self.__filter_blast()

            total_hits = self.blast_df.groupby("sseqid").size().reset_index(name="total_hits")
            passed_hits = (
                self.blast_df[(self.blast_df.passed_filters)].groupby("sseqid").size().reset_index(name="passed_hits")
            )
            df = pd.merge(total_hits, passed_hits, on="sseqid", how="outer").fillna(0)

            self.__summarise_groups(df)


def create_html_output(
    total_hits: int,
    passed_hits: int,
    yaml_details: aletheiaseq_validate_yml.SchemaDetails,
    methods_dict: dict,
    groups_matched: list[str],
    loci_results: dict[str, SchemaGroupSummary],
    sample_id: str,
    output_directory: str,
):
    html_summary_folder = Path(output_directory) / "html_summary_files"
    if html_summary_folder.is_dir():
        warnings.warn(
            "HTML output directory already exists. Any outputs files already created will be overwritten where the file prefix matches.",
            UserWarning,
            stacklevel=2,
        )
    else:
        html_summary_folder.mkdir(parents=False)

    species_list = ", ".join(methods_dict["schema_group_definitions"].keys())

    blast_filters = ", ".join(methods_dict["blast_filters"]["filters"])

    matched_group_list = ", ".join(groups_matched) if len(groups_matched) > 0 else "no species"

    blast_summary = (
        f"The {yaml_details.schema_type} loci definition was met for {matched_group_list}.\n\n"
        f"{sample_id} was processed using {yaml_details.name} ({yaml_details.schema_type}) version {yaml_details.version}.\n"
        f"The schema contains the following groups: {species_list}.\n"
        f"Blast hits were filtered by: {blast_filters} and at least {methods_dict['blast_filters']['minimum_blast_hits']} hits were required.\n"
        f"A total of {total_hits} hits were in the blast output and {passed_hits} met the filter thresholds.\n"
    )

    blast_summary_fp = html_summary_folder / "blast_headline.txt"
    blast_summary_fp.write_text(blast_summary)

    for name, summary in loci_results.items():
        group_result_fp = html_summary_folder / f"{name.lower().replace(' ', '_')}_locus_results.csv"
        pd.DataFrame(summary.loci).sort_values(["required", "passed_hits"], ascending=[False, True]).to_csv(
            group_result_fp, index=0
        )


def create_analysis_fields(
    sample_id: str,
    yaml_details: aletheiaseq_validate_yml.SchemaDetails,
    yaml_methods_dict: dict,
    headline_result: str,
    loci_results: dict[str, SchemaGroupSummary],
    server: Literal["mscape", "synthscape", "devscape"],
) -> tuple[oa.OnyxAnalysis, int]:
    """Set up fields dictionary used to populate analysis table for onyx.
    Arguments:
        sample_id -- Climb ID for sample
        yaml_details -- schema details from YAML file used to process input data
        yaml_methods_dict -- filters and thresholds used to process blast results
        headline_result -- Short description of main result
        loci_results -- Dictionary containing results per schema group defined in YAML
        server -- Server code is running on, one of "mscape", "devscape" or "synthscape"
    Returns:
        onyx_analysis -- Class containing required fields for input to onyx
                         analysis table
        exitcode -- Exit code for checks - will be 0 if all checks passed, 1 if any checks failed
    """
    onyx_analysis = oa.OnyxAnalysis()
    onyx_analysis.add_analysis_details(
        analysis_name=f"ukhsa-aletheiaseq-{yaml_details.name.lower().replace(' ', '-')}",
        analysis_description=f"Process blast output for {yaml_details.schema_type} using {yaml_details.name} version {yaml_details.version}",
    )
    onyx_analysis.add_package_metadata(package_name="AletheiaSeq")

    tool_versions = {
        "yaml_version": yaml_details.version,
    }
    methods_versions_fail = onyx_analysis.add_versions_to_methods(tool_versions=tool_versions)

    methods_fail = onyx_analysis.add_methods(methods_dict=yaml_methods_dict)

    results_fail = onyx_analysis.add_results(
        top_result=headline_result, results_dict={name: asdict(summary) for name, summary in loci_results.items()}
    )
    onyx_analysis.add_server_records(sample_id=sample_id, server_name=server)
    required_field_fail, attribute_fail = onyx_analysis.check_analysis_object(publish_analysis=False)

    exitcode = 1 if any([methods_versions_fail, methods_fail, results_fail, required_field_fail, attribute_fail]) else 0

    return onyx_analysis, exitcode


def parse(parsed_args: dict):
    asp = AletheiaSeqParser(parsed_args.yaml)
    asp.process_blast(parsed_args.blast_out)

    total_hits = asp.blast_df.shape[0]
    passed_hits = asp.blast_df[(asp.blast_df["passed_filters"])].shape[0]

    methods_dict = asp.yml.summarise_methods()

    groups_matched = []
    headlines = []
    for name, summary in asp.loci_results.items():
        if summary.all_req and summary.min_loci_pass:
            groups_matched.append(name)
        headlines.append(summary.headline(name))

    create_html_output(
        total_hits,
        passed_hits,
        asp.yml.SchemaDetails,
        methods_dict,
        groups_matched,
        asp.loci_results,
        parsed_args.sample,
        parsed_args.out_folder,
    )

    if parsed_args.onyx_server:
        onyx_output, code = create_analysis_fields(
            parsed_args.sample,
            asp.yml.SchemaDetails,
            methods_dict,
            ";".join(headlines),
            asp.loci_results,
            parsed_args.onyx_server,
        )

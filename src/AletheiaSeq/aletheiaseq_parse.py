import logging

import pandas as pd

from AletheiaSeq import aletheiaseq_validate_yml


class AletheiaSeqParser:
    def __init__(self, yaml_file: str, out_folder: str, sample_id: str, onyx: bool):
        self.yml = aletheiaseq_validate_yml.validate_yaml(yaml_file)
        self.out_folder = out_folder
        self.sample_id = sample_id
        self.onyx = onyx
        self.loci_results = {}

    def __load_blast(self, blast_out: str):
        columns = self.yml.BlastDetails.outfmt
        self.blast_df = pd.read_table(blast_out, sep="\t", header=None, names=columns)
        logging.info(
            f"Loaded blast output file with column headers: {','.join(columns)}. Number of rows loaded: {self.blast_df.shape[0]}."
        )

    def __filter_blast(self):
        filters = self.yml.BlastDetails.filters

        if len(filters) > 0:
            self.blast_df["passed_filters"] = self.blast_df.eval(" & ".join(filters))
        else:
            logging.info("No filters are present in the YAML file so all rows are considered QC pass")
            self.blast_df["passed_filters"] = True

    def __report_negative(self):
        for group in self.yml.SchemaGroups:
            defining_loci_len = len(group.defining_loci)
            df = pd.DataFrame(
                {
                    "sseqid": group.defining_loci,
                    "total_hits": [0] * defining_loci_len,
                    "passed_hits": [0] * defining_loci_len,
                    "above_min_threshold": [False] * defining_loci_len,
                }
            )
            df["required"] = df.sseqid.isin(group.required_loci)

            all_req = False
            min_loci = 0
            min_loci_pass = False

            self.loci_results[group.name] = {
                "results_df": df,
                "all_req": all_req,
                "min_loci": min_loci,
                "min_loci_pass": min_loci_pass,
            }

    def __summarise_groups(self):
        for group in self.yml.SchemaGroups:
            df = self.blast_summary[(self.blast_summary.sseqid.isin(group.defining_loci))].copy()
            df = df.merge(pd.DataFrame({"sseqid": group.defining_loci}), on="sseqid", how="outer").fillna(0)
            df["above_min_threshold"] = df.passed_hits >= group.minimum_loci_required
            df["required"] = df.sseqid.isin(group.required_loci)

            all_req = set(group.required_loci).issubset(
                set(df[df.above_min_threshold].sseqid)
            )  # this should still evaluate as true even if there are no required loci
            min_loci = df[df.above_min_threshold].shape[0]
            min_loci_pass = min_loci >= group.minimum_loci_required

            self.loci_results[group.name] = {
                "results_df": df,
                "all_req": all_req,
                "min_loci": min_loci,
                "min_loci_pass": min_loci_pass,
            }

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
            self.blast_summary = pd.merge(total_hits, passed_hits, on="sseqid", how="outer").fillna(0)

            self.__summarise_groups()


def parse(parsed_args):
    asp = AletheiaSeqParser(parsed_args.yaml, parsed_args.out_folder, parsed_args.sample, parsed_args.onyx)
    asp.process_blast(parsed_args.blast_out)

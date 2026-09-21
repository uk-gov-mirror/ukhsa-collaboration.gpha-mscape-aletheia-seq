import argparse
import pathlib
import textwrap
import warnings

from AletheiaSeq import aletheiaseq_parse, aletheiaseq_prepare


def _is_valid_path(path: str) -> str:
    path_obj = pathlib.Path(path)

    if path_obj.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid path")


def _is_valid_blast(base_path: str) -> str:
    required_files = [".nhr", ".nin", ".nsq", ".njs"]

    for ftype in required_files:
        path_obj = pathlib.Path(base_path + ftype)

        if path_obj.is_file():
            continue
        else:
            raise argparse.ArgumentTypeError(f"BlastDB:{ftype} is not present in directory given")

    return base_path


def _check_folder(folder: str) -> str:
    path_obj = pathlib.Path(folder)

    if path_obj.is_dir():
        warnings.warn(
            "Output directory already exists. Any outputs files already created will be overwritten where the file prefix matches.",
            UserWarning,
            stacklevel=2,
        )
        return folder
    elif path_obj.is_file():
        raise argparse.ArgumentTypeError(f"outfolder:{folder} is a path not a directory")
    else:
        try:
            path_obj.mkdir(parents=False)
            return folder
        except FileNotFoundError:
            raise argparse.ArgumentTypeError(f"outfolder:{folder} is missing required parent folder(s).") from None
        except OSError:
            raise argparse.ArgumentTypeError(f"outfolder:{folder} is invalid or a system I/O error occurred.") from None


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="AletheiaSeq",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
                                       In Greek mythology, Aletheia is the divine personification and spirit of truth, sincerity, and disclosure. Her name translates literally to 'unconcealedness' or 'state of being unhidden'.
                                       This package is intended to replicate loci-based presence/absence confirmation used in reference laboratories for speciation and characterisation of species of interest as a way of producing a 'ground-truth' from sequence data.
                                    """),
    )

    # all methods require YAML file so add arg at top level
    parser.add_argument(
        "--yaml",
        "-y",
        dest="yaml",
        required=True,
        type=_is_valid_path,
        help="Path to YAML file containing loci and QC information required to run AletheiaSeq. See full documentation for YAMl file format.",
    )
    parser.add_argument(
        "--outfolder",
        "-o",
        dest="out_folder",
        required=True,
        type=_check_folder,
        help="Folder to store output files. Folder will be created if it does not exist, files will be overwritten if it does.",
    )

    # add subparsers for prepare and parse
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    prepare = subparsers.add_parser(
        "prepare",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=textwrap.dedent("""
                                                            Use the prepare subcommand to check your input files are in the correct format and to produce an example bash script for running blast and (optionally) skope.
                                                            AletheiaSeq will not run these commandline programs for you and you do not need to run this validation every time you want to use the programme.
                                                         """),
    )
    prepare.set_defaults(func=aletheiaseq_prepare.prepare)
    prepare.add_argument(
        "--fasta",
        "-f",
        dest="ref_fasta",
        required=True,
        type=_is_valid_path,
        help="Path to FASTA file containing reference sequences for loci of interest. N.B. This must be the fasta file used to make the blast database.",
    )
    prepare.add_argument(
        "--blast_db",
        "-db",
        dest="blast_db",
        required=True,
        type=_is_valid_blast,
        help="Path to blast db (filename prefix only).",
    )
    prepare.add_argument(
        "--skope_idx",
        "-s",
        dest="skope_idx",
        required=False,
        type=_is_valid_path,
        help="Path to Skope index. Index will not be validated but if provided the skope command will be included in the bash script generated.",
    )

    parse = subparsers.add_parser(
        "parse",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=textwrap.dedent("""
                                                          Use the parse subcommand to process blast and (optionally) skope outputs and produce a report summarising the loci detected against the criteria set out in the input yaml.
                                                          To produce a json compatible with the onyx analysis table you will need to use the --onyx flag.
                                                       """),
    )
    parse.set_defaults(func=aletheiaseq_parse.parse)
    parse.add_argument(
        "--sample_id",
        "-id",
        dest="sample",
        type=str,
        required=True,
        help="Relevant ID for sample being processed. If outputs are being ingested into onyx this should be the CLIMB ID.",
    )
    parse.add_argument(
        "--blast_out",
        "-b",
        dest="blast_out",
        required=True,
        type=_is_valid_path,
        help="Path to blast output. N.B. outfmt must match what is detailed in the YAML to parse the output correctly.",
    )
    parse.add_argument("--skope_out", required=False, type=_is_valid_path, help="Path to output of skope classify.")
    parse.add_argument(
        "--onyx_server",
        required=False,
        type=str,
        choices=["mscape", "synthscape", "devscape"],
        help="Onyx server for analysis record. If not provided no analysis table result will be generated.",
    )
    parse.add_argument(
        "--publish",
        "-p",
        dest="publish",
        required=False,
        action="store_true",
        default=False,
        help="Flag to indicate whether onyx analysis object should be pushed to the database",
    )
    # this stores a bool counterintuitively
    # the default is to store True so that all runs are uploaded with dryrun=True
    # unless the --no_dryrun flag is explicitly set, in which case it is dryrun=False
    parse.add_argument(
        "--no_dryrun",
        dest="dryrun",
        required=False,
        action="store_false",
        default=True,
        help="Flag to indicate that onyx upload should not be a test upload (default is to perform test upload to prevent accidental publishing)",
    )

    return parser

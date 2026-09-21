import logging

from Bio.SeqIO import parse

from AletheiaSeq import aletheiaseq_validate_yml


def get_seqs_from_fasta(ref_fasta: str) -> list[str]:
    ids = []
    for rec in parse(ref_fasta, "fasta"):
        ids.append(rec.id)

    return ids


def prepare(parsed_args):
    logging.info("Validating yaml file against required schema...")
    validated_yaml = aletheiaseq_validate_yml.validate_yaml(parsed_args.yaml)

    yaml_seq_list = set(validated_yaml.BlastDetails.seqs)

    fasta_seq_list = get_seqs_from_fasta(parsed_args.ref_fasta)

    if len(set(fasta_seq_list) - yaml_seq_list) != 0:
        message = "Additional sequences in reference fasta that are not defined in yaml"
        logging.error(message)
        raise ValueError(message)
    elif len(yaml_seq_list - set(fasta_seq_list)) != 0:
        message = "Additional sequences defined in yaml that are not present in reference fasta"
        logging.error(message)
        raise ValueError(message)
    else:
        logging.info("Sequences defined in yaml match sequences in reference fasta")

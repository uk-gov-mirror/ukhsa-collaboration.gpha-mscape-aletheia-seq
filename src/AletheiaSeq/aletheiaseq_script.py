import logging
import sys
from datetime import datetime
from pathlib import Path

from AletheiaSeq.aletheiaseq_commandline import argument_parser


def set_up_logger(log_filepath):
    """
    Set up logger, set to append mode so logs from older runs are not overwritten.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s")

    out_handler = logging.FileHandler(log_filepath, mode="a")
    out_handler.setFormatter(formatter)
    logger.addHandler(out_handler)

    return logger


def main():
    args = argument_parser().parse_args()

    today = datetime.today().strftime("%Y-%m-%d")
    set_up_logger(Path(args.out_folder) / f"aletheiaseq_{args.subcommand}_{today}.logfile")

    args.func(args)


if __name__ == "__main__":
    sys.exit(main())

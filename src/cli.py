import argparse
import sys

# TODO
#   add base path
#   add option for blast (local / net)
#   translate


def parse_cli():
    """
    :return: Возвращает словарь с ключами input_file_name, output_file_name и list_genes.
            input_file_name - входной файл
            output_file_name - файл на запись
            list_genes - список генов
    """
    input_file_name = None
    output_file_name = None
    list_genes = []

    parser = argparse.ArgumentParser(
        description="BLAST Picha pastoris genes (PAS_...) vs Saccharomyces cerevisiae"
    )

    parser.add_argument("args", nargs=argparse.REMAINDER)

    parser.add_argument(
        "stdin", nargs="?", type=argparse.FileType("r"), default=sys.stdin
    )

    parser.add_argument(
        "-f",
        "--file",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help='File with genes. Should be column "id"',
    )

    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output with blast result",
    )

    args = parser.parse_args().args
    stdin = parser.parse_args().stdin
    input_file = parser.parse_args().file
    output_file = parser.parse_args().output

    if input_file.name == "<stdin>":
        # Нету файла на вход
        if not sys.stdin.isatty():
            stdin = stdin.read().splitlines()
        else:
            stdin = []

        list_genes = args + stdin
    else:
        input_file_name = input_file.name

    if output_file.name != "<stdout>":
        output_file_name = output_file.name

    return {
        "input_file_name": input_file_name,
        "output_file_name": output_file_name,
        "list_genes": list_genes,
    }

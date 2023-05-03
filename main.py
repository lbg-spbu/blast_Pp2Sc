from src.cli import parse_cli
from src.parse_input import ParseInput
from src.pp2sc_blast import BlastPp2Sc
from src.sgd import add_sgd_data
from src.wright_output import GenerateOutput

# TODO
#   logging
#   save data into sqlite and restart (save state)
#   starting info

def main():
    parse_res = parse_cli()
    input_file = parse_res["input_file_name"]
    output_file = parse_res["output_file_name"]
    pichia_genes = parse_res["list_genes"]

    print("============ Start ============")

    if len(pichia_genes) == 0:
        pichia_genes = ParseInput(input_file=input_file).run()
    else:
        if any((not gene.startswith("PAS") for gene in pichia_genes)):
            raise ValueError("Гены должны начинаться на PAS")

    blast_result = BlastPp2Sc(pichia_genes=pichia_genes).run()

    whole_result = add_sgd_data(blast_result)

    GenerateOutput(result=whole_result, output=output_file).run()

    print("============= Done ============")


if __name__ == "__main__":
    main()

from typing import List

from bioservices import UniProt

from src.types import PichiaUniProt
from src.utils import retries


def get_uniprot_data(genes: List[str]) -> List[PichiaUniProt]:
    uniprot = UniProt()
    pp_uniprot_data: List[PichiaUniProt] = []

    @retries(times=6, delay=3)
    def _uniprot_search(query: str) -> str:
        return uniprot.search(
            query,
            frmt="tsv",
            columns="gene_names, protein_name, go_p, go_f, go_c",
        )

    for gene in genes:
        res_uniprot = _uniprot_search(gene)
        parsed_result = res_uniprot.split("\n")[1].split("\t")
        try:
            pp_uniprot_res = PichiaUniProt(
                pp_gene_name=gene,
                protein_name=parsed_result[1],
                go_bio_process=parsed_result[2],
                go_mol_function=parsed_result[3],
                go_cell_component=parsed_result[4],
            )
        except IndexError:
            pp_uniprot_res = PichiaUniProt(
                pp_gene_name=gene,
                protein_name=res_uniprot,
                go_bio_process="",
                go_mol_function="",
                go_cell_component="",
            )

        pp_uniprot_data.append(pp_uniprot_res)

    return pp_uniprot_data

from typing import List

from pygenome import sg

from src.types import PichiaBlastFullData, PichiaSGD
from src.utils import parse_gene_name


def get_sgd_data(pichia_blast_data: List[PichiaBlastFullData]) -> List[PichiaSGD]:
    whole_pichia_data = []

    for data in pichia_blast_data:
        sc_gene_name = parse_gene_name(data.sc_description)

        if sc_gene_name is None:
            sc_gene_name = data.sc_description

        try:
            sgd_res = sg.stdgene[sc_gene_name]
        except Exception:
            sgd_description = None
        else:
            sgd_description = sgd_res.short_description

        whole_pichia_data.append(
            PichiaSGD(
                pp_gene_name=data.pp_gene_name,
                sgd_description=sgd_description,
            )
        )

    return whole_pichia_data

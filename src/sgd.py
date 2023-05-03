import re
from typing import List, Optional

from pygenome import sg

from src.types import PichiaBlastFullData, PichiaFullData


def parse_gene_name(description: str) -> Optional[str]:
    pattern = r"(?P<gene_name>[a-zA-Z]{3}\d{1,2})"
    res = re.search(pattern=pattern, string=description)

    if res is None:
        return None

    return res.group()


def add_sgd_data(pichia_blast_data: List[PichiaBlastFullData]) -> List[PichiaFullData]:
    whole_pichia_data = []

    for data in pichia_blast_data:
        sc_gene_name = parse_gene_name(data.sc_description)

        if sc_gene_name is None:
            sc_gene_name = data.sc_description

        try:
            sgd_res = sg.stdgene[sc_gene_name]
        except KeyError:
            sgd_description = None
        else:
            sgd_description = sgd_res.short_description

        whole_pichia_data.append(
            PichiaFullData(
                pp_gene_name=data.pp_gene_name,
                pp_gene_id=data.pp_gene_id,
                pp_description=data.pp_description,
                pp_protein_id=data.pp_protein_id,
                sc_accession=data.sc_accession,
                sc_description=data.sc_description,
                query_cover=data.query_cover,
                e_value=data.e_value,
                per_ident=data.per_ident,
                acc_le=data.acc_le,
                sgd_description=sgd_description,
            )
        )

    return whole_pichia_data

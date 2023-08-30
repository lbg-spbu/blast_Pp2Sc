from dataclasses import dataclass
from typing import Optional


@dataclass
class PichiaXP:
    pp_gene_name: str
    pp_gene_id: str
    pp_description: str
    pp_protein_id: str


@dataclass
class SaccharomycesBlast:
    sc_accession: Optional[str]
    sc_description: Optional[str]
    query_cover: Optional[str]
    e_value: Optional[str]
    per_ident: Optional[float]
    acc_le: Optional[float]


@dataclass
class PichiaBlastFullData(PichiaXP, SaccharomycesBlast):
    pass


@dataclass
class PichiaSGD:
    pp_gene_name: str
    sgd_description: Optional[str]

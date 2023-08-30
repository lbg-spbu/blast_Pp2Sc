from typing import List, Optional

import pandas as pd
from tqdm import tqdm

from src.types import PichiaBlastFullData, PichiaSGD, PichiaUniProt


class GenerateOutput:
    def __init__(
        self,
        blast_result: List[PichiaBlastFullData],
        sgd_result: List[PichiaSGD],
        uniprot_result: List[PichiaUniProt],
        output: Optional[str] = None,
    ):
        self.blast_result = blast_result
        self.sgd_result = sgd_result
        self.uniprot_result = uniprot_result
        self.output = output

    def run(self):
        if (
            len(self.blast_result) != len(self.sgd_result)
            or len(self.blast_result) != len(self.uniprot_result)
            or len(self.uniprot_result) != len(self.sgd_result)
        ):
            raise ValueError(
                f"Length of results don't match. ({self.blast_result=},{self.sgd_result=},{self.uniprot_result=})"
            )

        if self.output is None:
            self._print_output()
            return

        if self.output.endswith("csv"):
            self._csv_output()
        elif self.output.endswith("tsv"):
            self._tsv_output()
        else:
            self._print_output()

    @staticmethod
    def _check_result_ids(blast_res: str, sgd_res: str, uniprot_res: str):
        if blast_res != sgd_res or blast_res != uniprot_res or uniprot_res != sgd_res:
            raise ValueError("Gene names do not match")

    def _pandas_save(self) -> pd.DataFrame:
        data = pd.DataFrame(
            {
                "Pichia gene name": [],
                "Pichia gene ID": [],
                "Pichia description": [],
                "Pichia protein ID": [],
                "Pichia protein name": [],
                "Pichia GO bio process": [],
                "Pichia GO mol function": [],
                "Pichia GO cell component": [],
                "Saccharomyces accession ID": [],
                "Saccharomyces protein name": [],
                "Query Cover": [],
                "E-value": [],
                "Per. Ident": [],
                "Acc. Len": [],
                "SGD description": [],
            }
        )

        for blast_res, sgd_res, uniprot_res in tqdm(
            zip(self.blast_result, self.sgd_result, self.uniprot_result), total=len(self.blast_result)
        ):
            self._check_result_ids(blast_res.pp_gene_name, sgd_res.pp_gene_name, uniprot_res.pp_gene_name)
            row = pd.DataFrame(
                {
                    "Pichia gene name": [blast_res.pp_gene_name],
                    "Pichia gene ID": [blast_res.pp_gene_id],
                    "Pichia description": [blast_res.pp_description],
                    "Pichia protein ID": [blast_res.pp_protein_id],
                    "Pichia protein name": [uniprot_res.protein_name],
                    "Pichia GO bio process": [uniprot_res.go_bio_process],
                    "Pichia GO mol function": [uniprot_res.go_mol_function],
                    "Pichia GO cell component": [uniprot_res.go_cell_component],
                    "Saccharomyces accession ID": [blast_res.sc_accession],
                    "Saccharomyces protein name": [blast_res.sc_description],
                    "Query Cover": [blast_res.query_cover],
                    "E-value": [blast_res.e_value],
                    "Per. Ident": [blast_res.per_ident],
                    "Acc. Len": [blast_res.acc_le],
                    "SGD description": [sgd_res.sgd_description],
                }
            )
            data = pd.concat([data, row], ignore_index=True, axis=0)
        return data

    def _print_output(self):
        print("Result:")
        print(
            "Pp_gene_name",
            "Pp_gene_ID",
            "Pp_description",
            "Pp_protein_ID",
            "Pp_protein_name",
            "Pp_GO_bio_process",
            "Pp_GO_mol_function",
            "Pp_GO_cell_component",
            "Sc_accessionID",
            "Sc_protein_name",
            "Query_Cover",
            "E-value",
            "Per_Ident",
            "Acc_Len",
            "SGD_description",
            sep="\t",
        )
        for blast_res, sgd_res, uniprot_res in zip(self.blast_result, self.sgd_result, self.uniprot_result):
            self._check_result_ids(blast_res.pp_gene_name, sgd_res.pp_gene_name, uniprot_res.pp_gene_name)
            print(
                blast_res.pp_gene_name,
                blast_res.pp_gene_id,
                blast_res.pp_description,
                blast_res.pp_protein_id,
                uniprot_res.protein_name,
                uniprot_res.go_bio_process,
                uniprot_res.go_mol_function,
                uniprot_res.go_cell_component,
                blast_res.sc_accession,
                blast_res.sc_description,
                blast_res.query_cover,
                blast_res.e_value,
                blast_res.per_ident,
                blast_res.acc_le,
                sgd_res.sgd_description,
                sep="\t",
            )

    def _csv_output(self):
        print("Generating output table in .csv")
        data = self._pandas_save()
        print("Saving table")
        data.to_csv(self.output, sep=",")

    def _tsv_output(self):
        print("Generating output table in .tsv")
        data = self._pandas_save()
        print("Saving table")
        data.to_csv(self.output, sep="\t")

from typing import List, Optional

import pandas as pd
from tqdm import tqdm

from src.types import PichiaBlastFullData, PichiaSGD


class GenerateOutput:
    def __init__(
        self,
        blast_result: List[PichiaBlastFullData],
        sgd_result: List[PichiaSGD],
        output: Optional[str] = None,
    ):
        self.blast_result = blast_result
        self.sgd_result = sgd_result
        self.output = output

    def run(self):
        if len(self.blast_result) != len(self.sgd_result):
            raise ValueError(f"Length of results don't match. ({self.blast_result=},{self.sgd_result=})")

        if self.output is None:
            self._print_output()
            return

        if self.output.endswith("csv"):
            self._csv_output()
        elif self.output.endswith("tsv"):
            self._tsv_output()
        else:
            self._print_output()

    def _pandas_save(self) -> pd.DataFrame:
        data = pd.DataFrame(
            {
                "Pichia gene name": [],
                "Pichia gene ID": [],
                "Pichia description": [],
                "Pichia protein ID": [],
                "Saccharomyces accession ID": [],
                "Saccharomyces protein name": [],
                "Query Cover": [],
                "E-value": [],
                "Per. Ident": [],
                "Acc. Len": [],
                "SGD description": [],
            }
        )

        for blast_result, sgd_result in tqdm(zip(self.blast_result, self.sgd_result), total=len(self.blast_result)):
            if blast_result.pp_gene_name != sgd_result.pp_gene_name:
                raise ValueError("Gene names do not match")

            row = pd.DataFrame(
                {
                    "Pichia gene name": [blast_result.pp_gene_name],
                    "Pichia gene ID": [blast_result.pp_gene_id],
                    "Pichia description": [blast_result.pp_description],
                    "Pichia protein ID": [blast_result.pp_protein_id],
                    "Saccharomyces accession ID": [blast_result.sc_accession],
                    "Saccharomyces protein name": [blast_result.sc_description],
                    "Query Cover": [blast_result.query_cover],
                    "E-value": [blast_result.e_value],
                    "Per. Ident": [blast_result.per_ident],
                    "Acc. Len": [blast_result.acc_le],
                    "SGD description": [sgd_result.sgd_description],
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
            "Sc_accessionID",
            "Sc_protein_name",
            "Query_Cover",
            "E-value",
            "Per_Ident",
            "Acc_Len",
            "SGD_description",
            sep="\t",
        )
        for blast_result, sgd_result in zip(self.blast_result, self.sgd_result):
            if blast_result.pp_gene_name != sgd_result.pp_gene_name:
                raise ValueError("Gene names do not match")
            print(
                blast_result.pp_gene_name,
                blast_result.pp_gene_id,
                blast_result.pp_description,
                blast_result.pp_protein_id,
                blast_result.sc_accession,
                blast_result.sc_description,
                blast_result.query_cover,
                blast_result.e_value,
                blast_result.per_ident,
                blast_result.acc_le,
                sgd_result.sgd_description,
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

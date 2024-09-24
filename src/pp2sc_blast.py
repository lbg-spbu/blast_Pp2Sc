import os
import tempfile
from io import StringIO
from time import sleep
from typing import Generator, List, Literal, Union

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

from settings import DB_BASE_FILE, EMAIL, TMP_FILES_FOLDER
from src.types import PichiaBlastFullData, PichiaXP, SaccharomycesBlast
from src.utils import retries


class BlastPp2Sc:
    # TODO description
    """
    1. Из PAS_... (P.P) получаем protein accession number (XP...)
    2. BLAST из 1. на S.cerevisiae
    3. Возвращаем данные - result
    """

    EMAIL = EMAIL

    def __init__(self, pichia_genes: List[str]):
        self.pichia_genes = pichia_genes
        self.full_result: List[PichiaBlastFullData] = []

    def run(self, option: Union[Literal["local"], Literal["net"]] = "local") -> List[PichiaBlastFullData]:
        try:
            self._protein_blast_pp2sc(option)
        finally:
            if option == "local":
                for tmp_file in os.listdir(TMP_FILES_FOLDER):
                    if tmp_file.strip() == ".gitkeep":
                        continue
                    os.remove(os.path.join(TMP_FILES_FOLDER, tmp_file))

        return self.full_result

    @staticmethod
    @retries(times=6, delay=3)
    def _entrez_esearch(term):
        with Entrez.esearch(db="gene", term=term) as handle:
            search_record = Entrez.read(handle)
        return search_record

    @staticmethod
    @retries(times=6, delay=3)
    def _entrez_efetch(db, id, retmode, **kwargs):
        with Entrez.efetch(db=db, id=id, retmode=retmode, **kwargs) as handle:
            gene_record = Entrez.read(handle)
        return gene_record

    def _get_xp_from_pas(self) -> Generator[PichiaXP, None, None]:
        # TODO description
        """ """
        if len(self.pichia_genes) == 0:
            raise ValueError("Список генов пустой")  # TODO translate

        Entrez.email = self.EMAIL

        for i, gene in enumerate(self.pichia_genes):
            # TODO logging
            print(f"{i + 1}/{len(self.pichia_genes)} - Process gene - '{gene}'")

            search_record = self._entrez_esearch(gene)
            gene_id = search_record["IdList"][0]
            gene_record = self._entrez_efetch(db="gene", id=gene_id, retmode="xml")

            try:
                protein_description = gene_record[0]["Entrezgene_prot"]["Prot-ref"]["Prot-ref_name"][0]
            except KeyError:
                try:
                    protein_description = gene_record[0]["Entrezgene_prot"]["Prot-ref"]["Prot-ref_desc"]
                except KeyError:
                    protein_description = None

            protein_info = gene_record[0]["Entrezgene_locus"][0]["Gene-commentary_products"][0][
                "Gene-commentary_products"
            ][0]

            assert protein_info["Gene-commentary_type"].attributes["value"] == "peptide"

            prot_accession = protein_info["Gene-commentary_accession"]

            yield PichiaXP(
                pp_gene_name=gene,
                pp_gene_id=gene_id,
                pp_description=protein_description,
                pp_protein_id=prot_accession,
            )

    @staticmethod
    def _parse_blast_xml(blast_result: StringIO) -> SaccharomycesBlast:
        # TODO description
        """ """
        for record in NCBIXML.parse(blast_result):
            if len(record.alignments) == 0:
                return SaccharomycesBlast(
                    sc_accession=None,
                    sc_description=None,
                    query_cover=None,
                    e_value=None,
                    per_ident=None,
                    acc_le=None,
                )

            filtered_alignments = list(filter(lambda x: "NP" in x.accession, (i for i in record.alignments)))

            if len(filtered_alignments) == 0:
                align = record.alignments[0]
            else:
                align = filtered_alignments[0]

            sc_accession = align.accession
            sc_description = align.hit_def
            accession_len = align.length

            for hsp in align.hsps:
                percent_identity = hsp.identities / hsp.align_length * 100
                query_cover = (hsp.query_end - hsp.query_start + 1) / record.query_length * 100
                e_value = hsp.expect

                return SaccharomycesBlast(
                    sc_accession=sc_accession,
                    sc_description=sc_description,
                    query_cover=query_cover,
                    e_value=e_value,
                    per_ident=percent_identity,
                    acc_le=accession_len,
                )

    @staticmethod
    def _get_pichia_protein_fasta(protein_id: str):
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        return str(record.seq)

    def _protein_blast_pp2sc(self, option: Union[Literal["local"], Literal["net"]]) -> None:
        for pichia_data in self._get_xp_from_pas():
            if option == "local":
                pichia_protein_fasta = self._get_pichia_protein_fasta(pichia_data.pp_protein_id)

                with tempfile.NamedTemporaryFile(mode="w", dir=TMP_FILES_FOLDER, delete=False) as tmp:
                    tmp.write(pichia_protein_fasta)
                    blastp_cline = NcbiblastpCommandline(
                        query=tmp.name,
                        db=DB_BASE_FILE,
                        outfmt=5,
                    )

                sleep(0.3)
                blast = StringIO(blastp_cline()[0])

            elif option == "net":
                sleep(0.3)
                blast = NCBIWWW.qblast(
                    program="blastp",
                    database="refseq_protein",
                    sequence=pichia_data.pp_protein_id,
                    entrez_query="txid559292[ORGN]",
                    expect=0.001,
                )
            else:
                raise ValueError("Option ...")  # TODO

            blast_result = self._parse_blast_xml(blast)

            self.full_result.append(
                PichiaBlastFullData(
                    pp_gene_name=pichia_data.pp_gene_name,
                    pp_gene_id=pichia_data.pp_gene_id,
                    pp_description=pichia_data.pp_description,
                    pp_protein_id=pichia_data.pp_protein_id,
                    sc_accession=blast_result.sc_accession,
                    sc_description=blast_result.sc_description,
                    query_cover=blast_result.query_cover,
                    e_value=blast_result.e_value,
                    per_ident=blast_result.per_ident,
                    acc_le=blast_result.acc_le,
                )
            )

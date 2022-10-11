from time import sleep

from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML


class BlastPp2Sc:
    """
    1. Из PAS_... (P.P) получаем protein accession number (XP...)
    2. BLAST из 1. на S.cerevisiae
    3. Возвращаем данные - result

    result = {
        "PAS_...": [
            PP_gene_id,
            PP_description,
            PP_protein_id,
            SC_Accession,
            SC_description,
            Query_Cover,
            E_value,
            Per_Ident,
            Acc_Len
        ],
        "PAS_...": [
            ...
        ]
    }
    """

    EMAIL = "antonsidorin@list.ru"

    def __init__(self, pichia_genes: list):
        self.pichia_genes = pichia_genes
        self.result = {}

    def get_xp_from_pas(self):
        if len(self.pichia_genes) == 0:
            raise ValueError("Список генов пустой")

        Entrez.email = self.EMAIL

        for i, gene in enumerate(self.pichia_genes):
            print(
                f"{i + 1}/{len(self.pichia_genes)} - Process xp_from_pas for gene - '{gene}'"
            )

            if self.result.get(gene) is not None:
                raise KeyError(f"Такой ген ({gene}) уже есть")

            self.result[gene] = []

            with Entrez.esearch(db="gene", term=gene) as handle:
                search_record = Entrez.read(handle)

            gene_id = search_record["IdList"][0]

            with Entrez.efetch(db="gene", id=gene_id, retmode="xml") as handle:
                gene_record = Entrez.read(handle)

            protein_name = gene_record[0]["Entrezgene_prot"]["Prot-ref"][
                "Prot-ref_name"
            ][0]
            protein_info = gene_record[0]["Entrezgene_locus"][0][
                "Gene-commentary_products"
            ][0]["Gene-commentary_products"][0]
            assert protein_info["Gene-commentary_type"].attributes["value"] == "peptide"

            prot_accession = protein_info["Gene-commentary_accession"]

            self.result[gene].append(gene_id)
            self.result[gene].append(protein_name)
            self.result[gene].append(prot_accession)

            yield gene, prot_accession

    @staticmethod
    def parse_blast_xml(blast_result):
        """
        :param blast_result:
        :return: (
                SC_Accession,
                SC_description,
                Query_Cover,
                E_value,
                Per_Ident,
                Acc_Len
            )
        """
        for record in NCBIXML.parse(blast_result):
            if len(record.alignments) == 0:
                return None

            filtered_alignments = list(
                filter(lambda x: "NP" in x.accession, (i for i in record.alignments))
            )

            if len(filtered_alignments) == 0:
                align = record.alignments[0]
            else:
                align = filtered_alignments[0]

            sc_accession = align.accession
            sc_description = align.hit_def
            accession_len = align.length

            for hsp in align.hsps:
                percent_identity = hsp.identities / hsp.align_length * 100
                query_cover = hsp.align_length / record.query_length * 100
                e_value = hsp.expect

                return (
                    sc_accession,
                    sc_description,
                    query_cover,
                    e_value,
                    percent_identity,
                    accession_len,
                )

    def protein_blast_pp2sc(self):
        for gene, pichia_xp in self.get_xp_from_pas():
            print(f"\tBlast gene {gene}")
            sleep(0.3)
            blast = NCBIWWW.qblast(
                program="blastp",
                database="refseq_protein",
                sequence=pichia_xp,
                entrez_query="txid559292[ORGN]",
                expect=0.001,
            )
            blast_result = self.parse_blast_xml(blast)

            if blast_result is None:
                self.result[gene].extend([None] * 6)
            else:
                self.result[gene].extend(blast_result)

    def run(self):
        self.protein_blast_pp2sc()
        return self.result

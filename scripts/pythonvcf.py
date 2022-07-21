#!/usr/bin/env python

import argparse
import gzip
import os
import sys

class Sneff_transcript:
    def __init__(self, ann):
        splitted_ann = ann.split("|")
        self.allele = splitted_ann[0]
        self.effect = splitted_ann[1]
        self.impact = splitted_ann[2]
        self.gene = splitted_ann[3]
        self.geneid = splitted_ann[4]
        self.feature = splitted_ann[5]
        self.featureid = splitted_ann[6]
        self.biotype = splitted_ann[7]
        self.rank = splitted_ann[8]
        # alias HGVS_DNA, CODON): Variant in HGVS (DNA) notation
        self.hgvs_c = splitted_ann[9]
        # (alias HGVS, HGVS_PROT, AA): Variant in HGVS (protein) notation
        self.hgvs_p = splitted_ann[10]
        self.cdna_pos = splitted_ann[11]
        self.cdna_len = splitted_ann[12]
        self.cds_pos = splitted_ann[13]
        self.cds_len = splitted_ann[14]
        self.aa_pos = splitted_ann[15]

    def __str__(self):
        return "Effect: {} Impact: {} Gene: {}".format(self.effect,self.impact,self.gene)

    def __repr__(self):
        return "Effect: {} Impact: {} Gene: {}".format(self.effect,self.impact,self.gene)


class Variant:
    def __init__(self, vcf_line):
        splitted_line = vcf_line.split()

        number_of_columns = len(splitted_line)
        # If there aren't enough columns
        if number_of_columns < 10:
            sys.stderr.write("The number of fields is {} which is less than 10\n".format(number_of_columns))
            sys.exit(2)

        self.chromosome = splitted_line[0]
        self.position = int(splitted_line[1])
        self.identifier = splitted_line[2]
        self.reference = splitted_line[3]
        self.alternative = splitted_line[4]
        # If there is a numerical quality or not
        try:
            self.quality = float(splitted_line[5])
        except:
            self.quality = splitted_line[5]
        self.filter = splitted_line[6]
        self.info = splitted_line[7]
        self.format = splitted_line[8]

        sample_positions = range(9,number_of_columns)
        self.samples = [splitted_line[i] for i in sample_positions]

        self.samples_stats = {}
        self.__build_samples_stats()

        self.info_dict = {}
        self.__build_info_dict()

        # ClinVar aggregates information about genomic variation and its relationship to human health.
        self.clnsig = None
        # FATHMM Functional analysis through hidden markov model HMM
        # D: Deleterious T: Tolerated lower values are more deleterious
        self.fathmm_pred = None
        # If a FATHMM_score is <=0.5 (or rankscore <=0.81415) the corresponding non-synonymous SNP is predicted as "D(AMAGING)"; otherwise it is predicted as "T(OLERATED)"
        self.fathmm_score = None
        # https://fathmm.biocompute.org.uk/fathmmMKL.htm
        # Predictions are given as p-values in the range [0, 1]: values above 0.5 are predicted 
        # to be deleterious, while those below 0.5 are predicted to be neutral or benign. 
        # P-values close to the extremes (0 or 1) are the highest-confidence predictions that yield the highest accuracy.
        self.fathmm_mkl_coding_score = None
        self.fathmm_mkl_coding_pred = None

        # MutationAssessor score predicts the functional impact of amino-acid substitutions
        # in proteins based on evolutionary conservation of the affected amino acid in protein 
        # homologs. The score ranges from -5.2 to 6.5.
        # N = neutral
        # L = low
        # M = medium
        # H = high
        self.mutationassessor_score = None
        self.mutationassessor_pred = None

        # Description : Combined Annotation Dependent Depletion (CADD) - 
        # tool for scoring the deleteriousness of single nucleotide variants as well as insertion/deletions variants in the human genome.
        # CADD predicts a continuous phred-like score that ranges from 1 to 99, higher values indicating more deleterious cases
        self.cadd_phred = None
        # DANN is a functional prediction score based on a deep neural network.
        # The score can range from 0 to 1, when higher values are more likely to be deleterious.
        self.dann_score = None
        # The PolyPhen-2 score predicts the possible impact of an amino acid substitution on 
        # the structure and function of a human protein. This score represents the probability that a substitution is damaging.
        self.polyphen2_hvar_score = None
        self.polyphen2_hvar_pred = None
        # A SIFT score predicts whether an amino acid substitution affects protein function.
        self.sift_score = None
        self.sift_pred = None
        # Predicted loss of function effects for this variant. Format: 
        # 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'"
        self.lof = None
        # Predicted nonsense mediated decay effects for this variant. Format: 
        # 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected
        self.nmd = None
        # Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities
        self.mqranksum = None

        self._get_variant_effects()

        self.gnomad_genome_all = None

        self._get_variant_gnomad_stats()

        try:
            self.snpeff_transcipt_list = self._get_snpeff_transcripts(self.info_dict.get("ANN"))
        except:
            self.snpeff_transcipt_list = []


    def __str__(self):
        return("Chromosome: {}\nPosition: {}".format(self.chromosome, self.position))

    def __repr__(self):
        return("Chromosome: {}\nPosition: {}".format(self.chromosome, self.position))

    def __build_samples_stats(self):
        format_keys = self.format.split(":")
        sample_number = 0
        for sample in self.samples:
            n = 0
            self.samples_stats[sample_number] = {}
            splitted_sample = sample.split(":")
            if len(format_keys) != len(splitted_sample):
                sys.exit("This format {} and this sample {} have a different number of variables".format(self.format, sample))
            for k in format_keys:
                self.samples_stats[sample_number][k] = splitted_sample[n]
                n += 1
            sample_number += 1

    def get_sample_stats(self):
        """
        Parameters
        ----------
        sample : self

        Returns
        -------
        str
            A tab separated values string with GT AD and DP values for every sample in the vcf
        """
        s = ""
        for sample in self.samples_stats.keys():
            if sample > 0:
                s = s + "\t"
            sample_data = self.samples_stats.get(sample)
            # The Genotype
            gt = sample_data.get("GT")
            if gt == None:
                gt = "."
            s = s + gt 
            # Allelic depths for the ref and alt alleles in the order listed
            ad = sample_data.get("AD")
            if ad == None:
                ad = "."
            s = s + "\t" + ad
            # Approximate read depth
            dp = sample_data.get("DP")
            if dp == None:
                dp = "."
            s = s + "\t" + dp
        return s

    def _get_snpeff_transcripts(self, ANN_field):
        snpeff_transcipt_list = []
        transcripts = ANN_field.split(",")
        for transcript in transcripts:
            snpeff_transcipt = Sneff_transcript(transcript)
            snpeff_transcipt_list.append(snpeff_transcipt)
        return snpeff_transcipt_list

    def get_sample_values(self, key_values_dictionary):
        d = {}
        # GT: Genotype
        gt = key_values_dictionary.get("GT")
        d["gt"] = gt
        # AD: Allelic depths for the ref and alt alleles in the ord er listed
        ad = key_values_dictionary.get("AD")
        if ad == None:
            ad = "."
        d["ad"] = ad
        # DP: Approximate read depth
        dp = key_values_dictionary.get("DP")
        if dp == None:
            dp = "."
        d["dp"] = dp
        frequence =  [int(freq) for freq in ad.split(",")]
        # If ad has only one variable use the RD value 
        # RD: Depth of reference-supporting bases
        if len(frequence) == 1:
            frequence = [int(ad), int(key_values_dictionary.get("RD"))]
        ref,alt = frequence[0],frequence[1]
        d["ref"] = ref
        d["alt"] = alt
        if frequence[0] != 0:
            tumoral_variant_frequency = str(round(float(alt) / (float(ref) + float(alt)), 5))
        else:
            tumoral_variant_frequency = "Unknown"
        d["tumoral_variant_frequency"] = tumoral_variant_frequency
        return d


    def get_transcript_table_lines(self, name = None):
        lines = []
        starting_line = "{}\t{}\t{}\t{}\t{}\t{}\t".format(self.chromosome,
                self.position,
                self.identifier,
                self.reference,
                self.alternative,
                self.filter)

        # A string consisting of the disease name used by the database specified by CLNDISDB
        if 'CLNDN' not in self.info_dict:
            clndn = '.'
        else:
            clndn = self.info_dict["CLNDN"]

        type_of_mutation = "."
        if "recessive" in clndn: type_of_mutation = "recessive"
        if "dominant" in clndn: type_of_mutation = "dominant"

        more_info = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                clndn, type_of_mutation, self.clnsig, self.gnomad_genome_all, self.fathmm_pred,
                self.fathmm_score, self.fathmm_mkl_coding_score, self.fathmm_mkl_coding_pred, self.cadd_phred, self.dann_score,
                self.polyphen2_hvar_score, self.polyphen2_hvar_pred, self.sift_score, self.sift_pred,
                self.lof, self.nmd, self.mqranksum)

        unnamed_columns = ""
        number_of_samples = len(self.samples_stats.keys())
        i = 0
        for sample in self.samples_stats.keys():
            values = self.samples_stats[sample]
            sample_stats = self.get_sample_values(values)
            unnamed_columns = unnamed_columns +  sample_stats.get("gt") + "\t" + sample_stats.get("ad") + "\t" + sample_stats.get("dp") + "\t" + str(sample_stats.get("ref")) + "\t" + str(sample_stats.get("alt")) + "\t" + sample_stats.get("tumoral_variant_frequency")
            i += 1
            if i <= number_of_samples:
                unnamed_columns = unnamed_columns + "\t"


        snpeff_transcipts = []
        if len(self.snpeff_transcipt_list) != 0:
            for transcript in self.snpeff_transcipt_list:
                snpeff_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(transcript.effect,
                        transcript.impact, transcript.gene, transcript.geneid,
                        transcript.biotype, transcript.hgvs_c, transcript.hgvs_p,
                        transcript.cdna_pos, transcript.cds_pos, transcript.aa_pos)
                snpeff_transcipts.append(snpeff_line)
        else:
            snpeff_line = ".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t"
            snpeff_transcipts.append(snpeff_line)
        for transcript in snpeff_transcipts:
            complete_line = starting_line + transcript + more_info + "\t" + unnamed_columns
            lines.append(complete_line)
        return lines

    def get_genotype_field(self, sample, field):
        """
        Parameters
        ----------
        sample : int
            The number of the sample
        field : the field to be returned (example DP)

        Returns
        -------
        str
            the asked field
        """
        try:
            self.samples_stats[sample - 1][field]
        except:
            return None
        return self.samples_stats[sample - 1][field]


    def __build_info_dict(self):
        splitted_info = self.info.split(";")
        for sub_info in splitted_info:
            splitted_sub_info = sub_info.split("=")
            if len(splitted_sub_info) == 2:
                self.info_dict[splitted_sub_info[0]] = splitted_sub_info[1]

    def _get_variant_effects(self):
        """
        Description 
        -----------
            Populates the Variant with prediction of effects
            Available prediction softwares: FATHMM, ClinVar (CLNSIG)
        """
        # FATHMM Functional analysis through hidden markov model HMM
        # D: Deleterious T: Tolerated lower values are more deleterious
        self.fathmm_pred = self.info_dict.get("FATHMM_pred")
        # If a FATHMM_score is <=-1.5 (or rankscore <=0.81415) the corresponding non-synonymous SNP is predicted as "D(AMAGING)"; otherwise it is predicted as "T(OLERATED)"
        self.fathmm_score = self.info_dict.get("FATHMM_score")
        # https://fathmm.biocompute.org.uk/fathmmMKL.htm
        # Predictions are given as p-values in the range [0, 1]: values above 0.5 are predicted 
        # to be deleterious, while those below 0.5 are predicted to be neutral or benign. 
        # P-values close to the extremes (0 or 1) are the highest-confidence predictions that yield the highest accuracy.
        self.fathmm_mkl_coding_score = self.info_dict.get("fathmm-MKL_coding_score")
        self.fathmm_mkl_coding_pred = self.info_dict.get("fathmm-MKL_coding_pred")
        # ClinVar aggregates information about genomic variation and its relationship to human health.
        # MutationAssessor score predicts the functional impact of amino-acid substitutions
        # in proteins based on evolutionary conservation of the affected amino acid in protein 
        # homologs. The score ranges from -5.2 to 6.5.
        # N = neutral
        # L = low
        # M = medium
        # H = high
        self.mutationassessor_score = self.info_dict.get("MutationAssessor_score")
        self.mutationassessor_pred = self.info_dict.get("MutationAssessor_pred")
        # CADD predicts a continuous phred-like score that ranges from 1 to 99, higher values indicating more deleterious cases
        self.cadd_phred = self.info_dict.get("CADD_phred")
        # DANN is a functional prediction score based on a deep neural network.
        # The score can range from 0 to 1, when higher values are more likely to be deleterious.
        self.dann_score = self.info_dict.get("DANN_score")
        # The PolyPhen-2 score predicts the possible impact of an amino acid substitution on 
        # the structure and function of a human protein. This score represents the probability that a substitution is damaging.
        # Ion Reporter Software reports the pph2-prob PolyPhen-2 score.
        # The PolyPhen-2 score ranges from 0.0 (tolerated) to 1.0 (deleterious). Variants with scores of 0.0 are predicted to be benign. Values closer to 1.0 are more confidently predicted to be deleterious. The score can be interpreted as follows:
        # -  0.0 to 0.15 -- Variants with scores in this range are predicted to be benign.
        # - 0.15 to 1.0 -- Variants with scores in this range are possibly damaging.
        # - 0.85 to 1.0 -- Variants with scores in this range are more confidently predicted to be damaging.
        # PolyPhen-2 and SIFT scores use the same range, 0.0 to 1.0, but with opposite meanings. A variant with a PolyPhen-2 score of 0.0 is predicted to be benign. A variant with a SIFT score of 1.0 is predicted to be benign.
        self.polyphen2_hvar_score = self.info_dict.get("Polyphen2_HDIV_score")
        self.polyphen2_hvar_pred = self.info_dict.get("Polyphen2_HVAR_pred")
        # A SIFT score predicts whether an amino acid substitution affects protein function.
        # The SIFT score ranges from 0.0 (deleterious) to 1.0 (tolerated). The score can be interpreted as follows:
        # - 0.0 to 0.05 -- Variants with scores in this range are considered deleterious. Variants with scores closer to 0.0 are more confidently predicted to be deleterious.
        # - 0.05 to 1.0-- Variants with scores in this range are predicted to be tolerated (benign). Variants with scores very close to 1.0 are more confidently predicted to be tolerated.
        self.sift_score = self.info_dict.get("SIFT_score")
        self.sift_pred = self.info_dict.get("SIFT_pred")
        # Predicted loss of function effects for this variant. Format:
        # 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'"
        self.lof = self.info_dict.get("LOF")
        if self.lof != None:
            self.lof = self.lof[1:len(self.lof)][:-1].split("|")[-1]
        # Predicted nonsense mediated decay effects for this variant. Format: 
        # 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected
        self.nmd = self.info_dict.get("NMD")
        if self.nmd != None:
            self.nmd = self.nmd[1:len(self.nmd)][:-1].split("|")[-1]
        self.mqranksum = self.info_dict.get("MQRankSum")

        self.clnsig = self.info_dict.get("CLNSIG")

    def _get_variant_gnomad_stats(self):
        if self.info_dict.get("gnomAD_genome_ALL", ".") == ".":
            self.gnomad_genome_all = 0
        else:
            self.gnomad_genome_all = float(self.info_dict.get("gnomAD_genome_ALL"))




def main():
    parser = argparse.ArgumentParser(description="Parse a VCF file")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf to be parsed", required=True)
    parser.add_argument('-n', '--name', action='store', type=str, help="The sample name", required=False, default=None)
    parser.add_argument('-f', '--first', action='store_true', help="Print only the first variant")
    parser.add_argument('-t', '--no_header', action='store_true', help="Do not print the header")
    args = parser.parse_args()

    vcf = args.vcf
    name = args.name
    first = args.first
    no_header = args.no_header

    vcf_exists = os.path.isfile(vcf)
    if vcf_exists == False:
        sys.stderr.write("File {} do not exists\n".format(vcf))
        sys.exit(2)


    trait = '.'

    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        header = "chromosome\tposition\tidentifier\treference\talternative\tfilter\t\
                effect\timpact\tgene\tgene_id\tbiotype\thgvs_c\thgvs_p\tcdna_pos\t\
                cds_pos\taa_pos\tCLNDN\ttype_of_mutation\tclnsig\tgnomad_freq\t\
                fathmm_pred\tfathmm_score\tfathmm_mkl_coding_score\tfathmm_mkl_coding_pred\t\
                cadd_phred\tdann_score\tpolyphen2_hvar_score\tpolyphen2_hvar_pred\t\
                sift_score\tsift_pred\tlof\tnmd\tmqranksum"
        if args.name:
            header = "sample_name" + "\t" + header
        if no_header != True:
            print(header)
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                variant = Variant(line)

                if 'LRT_pred' not in variant.info_dict:
                    LRT_pred = '.'
                else:
                    LRT_pred = variant.info_dict["LRT_pred"]
                sample0_gt = variant.samples_stats[0]["GT"]
                if sample0_gt == "0/1":
                    trait = "heterozygous"
                if sample0_gt == "1/1":
                    trait = "homozygous"
                variants = variant.get_transcript_table_lines(args.name)
                for v in variants:
                    if name:
                        v = name + "\t" + v
                    print(v)
                    if first == True:
                        break

if __name__ == "__main__":
    main()

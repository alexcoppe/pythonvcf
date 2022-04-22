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
        self.hgvs_c = splitted_ann[9]
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

        self._get_variant_effects()

        self.gnomad_genome_all = None

        self._get_variant_gnomad_stats()

        try:
            self.snpeff_transcipt_list = self._get_snpeff_transcripts(self.info_dict.get("ANN"))
            #print(self.snpeff_transcipt_list)
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

    def _get_snpeff_transcripts(self, ANN_field):
        snpeff_transcipt_list = []
        transcripts = ANN_field.split(",")
        for transcript in transcripts:
            snpeff_transcipt = Sneff_transcript(transcript)
            snpeff_transcipt_list.append(snpeff_transcipt)
        return snpeff_transcipt_list
        

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
        #print(self.fathmm_mkl_coding_pred)
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

        self.clnsig = self.info_dict.get("CLNSIG")

    def _get_variant_gnomad_stats(self):
        if self.info_dict.get("gnomAD_genome_ALL") == ".":
            self.gnomad_genome_all = 0
        else:
            self.gnomad_genome_all = float(self.info_dict.get("gnomAD_genome_ALL"))


def main():
    parser = argparse.ArgumentParser(description="Parse a VCF file")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf to be parsed", required=True)
    args = parser.parse_args()

    vcf = args.vcf

    vcf_exists = os.path.isfile(vcf)
    if vcf_exists == False:
        sys.stderr.write("File {} do not exists\n".format(vcf))
        sys.exit(2)

    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                variant = Variant(line)
                #print(variant.fathmm_score, variant.fathmm_pred, variant.fathmm_mkl_coding_score, variant.clnsig)
                #print(variant.gnomad_genome_all)
                #print(variant.info_dict.get("ANN"))
                #print("{} {}".format(variant.mutationassessor_pred, variant.mutationassessor_score))

                clndn = variant.info_dict["CLNDN"]
                LRT_pred = variant.info_dict["LRT_pred"]
                #print(LRT_pred)
                #print(variant.info_dict.get("Polyphen2_HDIV_score"))
                #print(variant.filter)
                sample0_gt = variant.samples_stats[0]["GT"]
                type_of_mutation = ""
                if sample0_gt == "0/1":
                    trait = "heterozygous"
                if sample0_gt == "1/1":
                    trait = "homozygous"
                if "recessive" in clndn: type_of_mutation = "recessive"
                if "dominant" in clndn: type_of_mutation = "dominant"
                if type_of_mutation != "":
                    to_print = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(trait, type_of_mutation, variant.chromosome, variant.position, variant.identifier, variant.reference, variant.alternative, sample0_gt, LRT_pred, clndn)
                    print(to_print)
                    #print("\n")

if __name__ == "__main__":
    main()

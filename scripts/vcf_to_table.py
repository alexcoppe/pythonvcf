#!/usr/bin/env python

import argparse
import gzip
import os
import sys


class Sneff_transcript:
    """
    A class used to represent a SnpEff transcript

    ...

    Attributes
    ----------
    allele : str
        Allele (or ALT)
    effect : str
        Annotation (a.k.a. effect): Annotated using Sequence Ontology terms 
    impact : str
       Putative_impact: A simple estimation of putative impact / deleteriousness : {HIGH, MODERATE, LOW, MODIFIER} 
    """ 
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
    """
    A class used to represent a Variant

    ...

    Attributes
    ----------
    chromosome : str
        the chromosome
    position : int
        the position of the variant
    identifier : str
        an identifier
    reference : str
        the bases on the refrence genome
    alternative : str
        the bases on the currently considered living thing
    quality : float
        Phred-scaled quality score for the assertion made in ALT
    filter : str
        filter status: PASS if this position has passed all filters
    info : str
        additional information
    info_dict : dictionary
        a dictionary containing all the key values from the info field

    """
    def __init__(self, vcf_line):
        """
        Parameters
        ----------
        vcf_line : string
            a line from a VCF file
        """

        self.splitted_line = vcf_line.split()

        number_of_columns = len(self.splitted_line)
        # If there aren't enough columns
        if number_of_columns < 8:
            sys.stderr.write("The number of fields is {} which is less than 10\n".format(number_of_columns))
            sys.exit(2)

        self.chromosome = self.splitted_line[0]
        self.position = int(self.splitted_line[1])
        self.identifier = self.splitted_line[2]
        self.reference = self.splitted_line[3]
        self.alternative = self.splitted_line[4]
        # If there is a numerical quality or not
        try:
            self.quality = float(self.splitted_line[5])
        except:
            self.quality = self.splitted_line[5]
        self.filter = self.splitted_line[6]
        self.info = self.splitted_line[7]

        self.info_dict = {}
        self.__build_info_dict()

        try:
            self.snpeff_transcipt_list = self.__get_snpeff_transcripts(self.info_dict.get("ANN"))
        except:
            self.snpeff_transcipt_list = []

        self.samples = []

        if number_of_columns > 8:
            self.__build_list_of_samples()

    def __build_list_of_samples(self):
        number_of_elements_from_split = len(self.splitted_line)
        format_elements = self.splitted_line[8].split(':')
        for i in range(9,number_of_elements_from_split):
            splitted_sample = self.splitted_line[i].split(':')
            if len(splitted_sample) != len(format_elements):
                print("Error")
            else:
                splitted_sample = self.splitted_line[i].split(':')
                d = {}
                i = 0
                for key in format_elements:
                    d[key] = splitted_sample[i]
                    i += 1
                self.samples.append(d)


    def __get_snpeff_transcripts(self, ANN_field):
        """
        Description 
        -----------
            A private method that get the SnpEff transcripts
        """
        snpeff_transcipt_list = []
        transcripts = ANN_field.split(",")
        for transcript in transcripts:
            snpeff_transcipt = Sneff_transcript(transcript)
            snpeff_transcipt_list.append(snpeff_transcipt)
        return snpeff_transcipt_list

    def __build_info_dict(self):
        """
        Description 
        -----------
            Populates the info_dict dictionary
        """
        splitted_info = self.info.split(";")
        for sub_info in splitted_info:
            splitted_sub_info = sub_info.split("=")
            if len(splitted_sub_info) == 2:
                if splitted_sub_info[1].isnumeric():
                    self.info_dict[splitted_sub_info[0]] = int(splitted_sub_info[1])
                elif splitted_sub_info[1].replace(".", "").isnumeric():
                    self.info_dict[splitted_sub_info[0]] = float(splitted_sub_info[1])
                else:
                    self.info_dict[splitted_sub_info[0]] = splitted_sub_info[1]



    def __str__(self):
        return("Chromosome: {}\nPosition: {}".format(self.chromosome, self.position))

    def __repr__(self):
        return("Chromosome: {}\nPosition: {}".format(self.chromosome, self.position))



def main():
    parser = argparse.ArgumentParser(description="Gets a VCF and shows a tsv with selected fields and subfields")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf to be parsed", required=True)
    parser.add_argument('-f', '--fields', action='store', type=str, help="The file containing the list of wanted fields and subfields", required=True)

    args = parser.parse_args()

    vcf = args.vcf
    wanted_fields_file_name = args.fields

    vcf_exists = os.path.isfile(vcf)
    if vcf_exists == False:
        sys.stderr.write("File {} do not exists\n".format(vcf))
        sys.exit(2)

    wanted_fields = []
    with (open(wanted_fields_file_name) as  wanted_fields_content):
        for line in wanted_fields_content:
            wanted_fields.append(line[:-1])

    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')

            if line[0] != '#':
                number_of_elements_in_line = len(line.split("\t"))
                if number_of_elements_in_line >= 8:
                    variant = Variant(line)
                    number_of_snpeff = len(variant.snpeff_transcipt_list)
                    s = ""
                    l = []
                    for transcript in range(0,number_of_snpeff):
                        l.append("")

                    for wanted_field in wanted_fields:
                        if ':' in wanted_field:
                            wf = wanted_field.split(':')[1]

                            if wf == "allele":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.allele + "\t"
                                    i += 1
                            if wf == "effect":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.effect + "\t"
                                    i += 1
                            if wf == "impact":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.impact + "\t"
                                    i += 1
                            if wf == "gene":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.gene + "\t"
                                    i += 1
                            if wf == "geneid":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.geneid + "\t"
                                    i += 1
                            if wf == "feature":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.feature + "\t"
                                    i += 1
                            if wf == "featureid":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.featureid + "\t"
                                    i += 1
                            if wf == "biotype":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.biotype + "\t"
                                    i += 1
                            if wf == "rank":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.rank + "\t"
                                    i += 1
                            if wf == "hgvs_c":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.hgvs_c + "\t"
                                    i += 1
                            if wf == "hgvs_p":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.hgvs_p + "\t"
                                    i += 1
                            if wf == "cdna_pos":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.cdna_pos + "\t"
                                    i += 1
                            if wf == "cdna_len":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.cdna_len + "\t"
                                    i += 1
                            if wf == "cds_pos":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.cds_pos + "\t"
                                    i += 1
                            if wf == "cds_len":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.cds_len + "\t"
                                    i += 1
                            if wf == "aa_pos":
                                i = 0
                                for transcript in variant.snpeff_transcipt_list:
                                    l[i] = l[i] + transcript.aa_pos + "\t"
                                    i += 1
                        elif ';' in wanted_field:
                            wf = wanted_field.split(';')[1]
                            for i in range(0,number_of_snpeff):
                                l[i] = l[i] + str(variant.info_dict[wf])  + "\t" 

                        elif '|' in wanted_field:
                            wf = wanted_field.split('|')[1]
                            s = ""
                            dff = 0
                            for dictionary in variant.samples:
                                value = dictionary.get(wf)
                                s = s + value  + "\t"
                            for i in range(0,number_of_snpeff):
                                l[i] = l[i] + s + "\t" 
                        else:
                            wf = wanted_field
                            if wf == "position":
                                for i in range(0,number_of_snpeff):
                                    l[i] = l[i] + str(variant.position) + "\t" 
                            elif wf == "chromosome":
                                for i in range(0,number_of_snpeff):
                                    l[i] = l[i] + str(variant.chromosome) + "\t" 
                            elif wf == "identifier":
                                for i in range(0,number_of_snpeff):
                                    l[i] = l[i] + str(variant.identifier) + "\t" 
                            elif wf == "reference":
                                for i in range(0,number_of_snpeff):
                                    l[i] = l[i] + str(variant.reference) + "\t" 
                            elif wf == "alternative":
                                for i in range(0,number_of_snpeff):
                                    l[i] = l[i] + str(variant.alternative) + "\t" 
                            elif wf == "quality":
                                for i in range(0,number_of_snpeff):
                                    l[i] = l[i] + str(variant.quality) + "\t" 
                            elif wf == "filter":
                                for i in range(0,number_of_snpeff):
                                    l[i] = l[i] + str(variant.filter) + "\t" 
                            elif wf == "info":
                                for i in range(0,number_of_snpeff):
                                    l[i] = l[i] + str(variant.info) + "\t" 


                    for k in l:
                        print(k)

                            
if __name__ == "__main__":
    main()

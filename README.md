# hdose

Harmonize dosage, not effect, for mendelian randomization.

harmonize-dosage is a program to harmonize the dosage and dosage info files to
the GWAS file in a Mendelian Randomization experiment. It accepts both dosage
files and vcf.gz/vcf/bcf files as input. In the latter case, a dosage file is
created from the vcf.

## Changelog

```
0.0.4 (30/10/11)
- much renaming and improving
0.0.3 (23/10/11)
- critical fix: used inverse AF for effect allele
```

## Install

```bash
pip install hdose
```

## Requirements

To use vcf/vcf.gz/bcf-files, bcftools is required on the path.

## Quick-start

```bash
# get the test data
git clone git@github.com:hunt-genes/harmonize_dosage.git hdose

# go to folder with test data
cd hdose

# view test data
tail -n 3 tests/data/*
# ==> tests/data/exposure.txt <==
# "#CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO" "FORMAT" "rs" "AF" "MAF" "RSQ" "GT"
# "1" Z 15 "Z:15_G/A" "G" "A" "." "PASS" "AF=0.285;MAF=0.36146;R2=0.99242" "GT:DS" "bs1" "0.63854" "0.36146" "0.99242" "IMPUTED"
# "2" Z 42 "Z:42_G/T" "G" "T" "." "PASS" "AF=0.8456;MAF=0.14664;R2=0.99933;ER2=0.97741;GENOTYPED" "GT:DS" "bs3" "0.14664" "0.14664" "0.99933" "GENOTYPED"
#
# ==> tests/data/outcome.txt <==
# SNP_hg18      SNP_hg19       rsid A1 A2   beta     se         N    P.value FreqA2
# chrZ:15 chrZ:15  bs1  A  g 0.0576 0.0071 123.02  1.982e-29           0.28500
# chrZ:42 chrZ:42 bs3  T  c 0.0191 0.0055 456.97  0.750e-33           0.84560
#
# ==> tests/data/tiny.vcf <==
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2	NA3
# Z	15	bs1	C	G	.	PASS	AF=0.123;MAF=0.456;R2=0.0;ER2=0.98798;GENOTYPED	GT:DS	1|1:1.300	1|1:0.010	0|1:2.000
# Z	42	bs3	T	G	.	PASS	AF=0.123;MAF=0.456;R2=0.0;ER2=0.99968;GENOTYPED	GT:DS	0|1:2.000	0|1:0.000	1|0:1.000

# run hdose
hdose -d tests/data/exposure.txt -v tests/data/tiny.vcf  -g tests/data/outcome.txt \
-o outfolder_quickstart -G 'A1,A2,rsid,FreqA1' -D 'REF,ALT,rs,AF'

# swapped 50.0
# incompatible 50.0

# view output files
# only one column in dosage file since incompatible snps are removed
head outfolder_quickstart/*
# ==> outfolder_quickstart/harmonized_dosage.csv <==
# ID bs1_dose
# NA1 0.7
# NA2 1.99
# NA3 0.0
#
# ==> outfolder_quickstart/harmonized_dosage_info.csv <==
# #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT rs AF MAF RSQ GT
# bs1 Z 15 Z:15_G/A A G . PASS AF=0.285;MAF=0.36146;R2=0.99242 GT:DS 0.63854 0.36146 0.9924200000000001 IMPUTED
#
# ==> outfolder_quickstart/missing_data.csv <==
#  rsid D1 D2 G1 G2 GAF DAF
#
# ==> outfolder_quickstart/snp_classifications.csv <==
# rsid D1 D2 G1 G2 GAF DAF SNPType
# bs1 G A A G 0.285 0.63854 swapped
# bs3 G T T C 0.8456 0.14664000000000002 incompatible
```

## Usage

(Might be slightly out of date.)

```
hdose
Harmonize the dosages between dosage_info and GWAS alleles for Mendelian Randomization.
(Visit github.com/endrebak/harmonize-dosage for examples and help.)
Usage:
    hdose --dosage-info=DSI --GWAS=GWS [--dosage=DOS] --output-dir=DIR
          [--tolerance=TOL] [--GWAS-columns=GWC] [--dosage-info-columns=DSC]
    hdose --help
Arguments:
    -d DSI --dosage-info DSI          Dosage info file
    -g GWS --GWAS GWS                 GWAS info file
    -v DOS --dosage DOS               Dosage file (bcf/dosage/vcf/vcf.gz)
    -t TOL --tolerance TOL            Tolerance for considering AFs similar.  [default: 0.08]
    -o DIR --output-dir DIR           Dir to put output files in.
    -D DSC --dosage-info-columns DSC  Comma-separated names of the columns containing the
                                      reference allele, alternate  allele, rsid and
                                      allele frequency in the dosage-info file. The allele frequency
                                      should belong to the alternate allele.
                                      [default: REF,ALT,rs,AF]
    -G GWC --GWAS-columns GWC         Comma-separated names of the columns containing the
                                      reference allele, alternate allele, rsid and
                                      allele frequency in the GWAS file. The allele frequency
                                      should belong to the alternate allele.
                                      [default: A1,A2,rsid,FreqA2]
Options:
    -h --help                 show this help message
```

## Contributing

If you find any bugs, please report them on the issues page.

## License

MIT

## Authors

Endre Bakken Stovner and Ben Michael Brumpton

## Thanks

* Humaira Rasheed (testing the script)

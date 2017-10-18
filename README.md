# hdose

Harmonize dosage, not effect, for mendelian randomization.

harmonize-dosage is a program to harmonize the dosage and exposure files to the
outcome file in a Mendelian Randomization experiment. It accepts both dosage and
vcf.gz/vcf/bcf files as input. In the latter case, a dosage file is created from
the vcf.

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
==> tests/data/exposure.txt <==
"#CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO" "FORMAT" "rs" "AF" "MAF" "RSQ" "GT"
"1" Z 15 "Z:15_G/A" "G" "A" "." "PASS" "AF=0.285;MAF=0.36146;R2=0.99242" "GT:DS" "bs1" "0.63854" "0.36146" "0.99242" "IMPUTED"
"2" Z 42 "Z:42_G/T" "G" "T" "." "PASS" "AF=0.8456;MAF=0.14664;R2=0.99933;ER2=0.97741;GENOTYPED" "GT:DS" "bs3" "0.14664" "0.14664" "0.99933" "GENOTYPED"

==> tests/data/outcome.txt <==
SNP_hg18      SNP_hg19       rsid A1 A2   beta     se         N    P.value FreqA1
chrZ:15 chrZ:15  bs1  A  g 0.0576 0.0071 123.02  1.982e-29           0.28500
chrZ:42 chrZ:42 bs3  T  c 0.0191 0.0055 456.97  0.750e-33           0.84560

==> tests/data/tiny.vcf <==
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2	NA3
Z	15	bs1	C	G	.	PASS	AF=0.123;MAF=0.456;R2=0.0;ER2=0.98798;GENOTYPED	GT:DS	1|1:1.300	1|1:0.010	0|1:2.000
Z	42	bs3	T	G	.	PASS	AF=0.123;MAF=0.456;R2=0.0;ER2=0.99968;GENOTYPED	GT:DS	0|1:2.000	0|1:0.000	1|0:1.000

# run hdose
hdose -e tests/data/exposure.txt -d tests/data/tiny.vcf  -o tests/data/outcome.txt \
-D outfolder_test -O 'A1,A2,rsid,FreqA1' -E 'REF,ALT,rs,AF'

# view output files
head outfolder_quickstart/*
==> outfolder_quickstart/harmonized_dosage.csv <==
ID bs1_dose
NA1 0.7
NA2 1.99
NA3 0.0

==> outfolder_quickstart/harmonized_exposure.csv <==
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT rs AF MAF RSQ GT
bs1 Z 15 Z:15_G/A A G . PASS AF=0.285;MAF=0.36146;R2=0.99242 GT:DS 0.63854 0.36146 0.9924200000000001 IMPUTED

==> outfolder_quickstart/missing_data.csv <==
 rsid E1 E2 O1 O2 OAF EAF

==> outfolder_quickstart/snp_classifications.csv <==
rsid E1 E2 O1 O2 OAF EAF SNPType
bs1 G A A G 0.285 0.63854 swapped
bs3 G T T C 0.8456 0.14664000000000002 incompatible
```

## Usage

```
hdose
Harmonize the dosages between exposure and outcome alleles for Mendelian Randomization.
(Visit github.com/endrebak/harmonize-dosage for examples and help.)
Usage:
    hdose --exposure=EXP --outcome=OUT --dosage=DOS --output-dir=DIR
          [--tolerance=TOL] [--outcome-columns=OC] [--exposure-columns=EC]
    hdose --help
Arguments:
    -e EXP --exposure EXP     Exposure file
    -o OUT --outcome OUT      Outcome file
    -d DOS --dosage DOS       Dosage file
    -t TOL --tolerance TOL    Tolerance for considering AFs similar.  [default: 0.08]
    -D DIR --output-dir DIR   Dir to put output files in.
    -E EC --exposure-columns  Comma-separated names of the columns containing the
                              reference allele, alternate allele, rsid and
                              allele frequency in the outcome file.
                              [default: REF,ALT,rs,AF]
    -O OC --outcome-columns   Comma-separated names of the columns containing the
                              reference allele, alternate allele, rsid and
                              allele frequency in the outcome file.
                              [default: A1,A2,rsid,FreqA1]
Options:
    -h --help                 show this help message
```

## Contributing

If you find any bugs, please report them on the issues page.

## License

MIT

## Authors

Endre Bakken Stovner and Ben Michael Brumpton

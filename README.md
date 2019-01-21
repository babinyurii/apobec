# apobec
### NGS data analysis for HBV crispr/cas9 research group 

## Requirements
- Python 3
- Biopython
- matplotlib
- seaborn
- pandas

## Description
This script takes fasta alignment as an input. The input file is the result of reads mapping onto the reference sequence. `count_snp_duplex.py` counts  SNP in dinucleotide duplex context.

`count_snp_duplex.py` outputs summary bar charts : 
![bars](output_example/bars.png)

and excel spreadsheets to further manipulate the data  :
![bars](output_example/raw_count_spread_sheet.PNG)

![bars](output_example/pivot_count.PNG)

![bars](output_example/pivot_percent.PNG)

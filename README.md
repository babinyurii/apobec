# apobec
`apobec` contains two in-house scripts that are used for NGS data analysis by HBV crispr/cas9 research group. 

## Installation
```
pip install apobec
```
## Usage
`apobec` can be used via shell or Jupyter Notebook. Create folder named `input_data` and put your fastas into it. Navigate into the directory which contains the `input_data` folder. Then import package :
```python
import apobec
```
and run via shell :
```
python -m apobec.count_snp_duplex
python -m apobec.snp_rate
```
or by Jupyter Notebook :
```python
%run -m apobec.count_snp_duplex
%run -m apobec.snp_rate
```

## Description
The scripts take fasta alignment as an input. The input file is the result of deep sequencing reads mapping onto the reference sequence and is imported from the Geneious software.

`count_snp_duplex.py` counts  SNP in dinucleotide duplex context.

`count_snp_duplex.py` outputs summary bar charts : 
![bars](output_example/bars.png)

and excel spreadsheets to further manipulate the data  :
![bars](output_example/raw_count_spread_sheet.PNG)

![bars](output_example/pivot_count.PNG)

![bars](output_example/pivot_percent.PNG)

`snp_rate.py` counts SNP in each read.

It outputs a distribution plot:
![dist](output_example/snp_dist.png)

and also raw count and summary statistics in excel spreadsheets:
![sum_stat](output_example/pivot_mutation_data.png)

## Requirements
- Python 3
- biopython
- matplotlib
- numpy
- pandas
- seaborn

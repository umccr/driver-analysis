# Driver analysis - Cancer Genome Interpreter

## Table of contents

<!-- vim-markdown-toc GFM -->
* [Description](#description)
* [Usage](#usage)
	* [Input data preparation](#input-data-preparation)
	* [Launching job](#launching-job)
	* [Downloading results](#downloading-results)
* [Output](#output)

<!-- vim-markdown-toc -->

<br>

### Description

[Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org) (CGI) is a tool able to flag validated oncogenic alterations and predict cancer drivers among mutations of unknown significance (see paper by [Tamborero et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/29592813)). Additionally, it flags genomic biomarkers of drug response with different levels of clinical relevance.


### Usage

[CGI](https://www.cancergenomeinterpreter.org) is freely available through a web interface at [http://www.cancergenomeinterpreter.org](http://www.cancergenomeinterpreter.org). The CGI resource can also be accessed programmatically by an **API** created via **REST** which is described below. Only registered users can make use of the API, since a **token** is needed for any communication between the end user and the REST API. It can be accessed in the following URL [https://www.cancergenomeinterpreter.org/api/v1](https://www.cancergenomeinterpreter.org/rest_api) and information about acceptable input format can be found at [https://www.cancergenomeinterpreter.org/formats](https://www.cancergenomeinterpreter.org/formats).

<br /> 

#### Input data preparation

[CGI](https://www.cancergenomeinterpreter.org/api/v1) REST API expects the input data to be in particular format, which is different from standard [MAF](https://software.broadinstitute.org/software/igv/MutationAnnotationFormat) file format. In order to **convert MAF into CGI**-comptaible files follow these steps

To run [CGI](https://www.cancergenomeinterpreter.org/api/v1) REST API one needs to ptovide the `job identifier`, `email` and `token`.

Here we assume that the input MAF file is `/data/example.maf`, e-mail address is `jacek.marzec@unimelb.edu.au`, the token is `abcd1234`, job ID is `123456789` and we refer to pancreatic cancer `PAAD`

```
# Reformat MAF file to be compatible with CGI
sed '/^#/ d' < /data/example.maf > /data/example.txt.tmp

3 Rename some column names
sed -i.bak 's/Chromosome/chr/g' /data/example.txt.tmp
sed -i.bak 's/Start_Position/pos/g' /data/example.txt.tmp
sed -i.bak 's/Strand/strand/g' /data/example.txt.tmp
sed -i.bak 's/Reference_Allele/ref/g' /data/example.txt.tmp
sed -i.bak 's/Tumor_Seq_Allele2/alt/g' /data/example.txt.tmp
sed -i.bak 's/Tumor_Sample_Barcode/sample/g' /data/example.txt.tmp
rm /data/example.txt.tmp.bak

## Keep only required columns to minimise the file size and speed up data transfer
awk -F'\t' -v cols=Hugo_Symbol,chr,pos,ref,alt,sample '(NR==1){n=split(cols,cs,",");for(c=1;c<=n;c++){for(i=1;i<=NF;i++)if($(i)==cs[c])ci[c]=i}}{for(i=1;i<=n;i++)printf "%s" FS,$(ci[i]);printf "\n"}' /data/example.txt.tmp > /data/example.txt

rm /data/example.txt.tmp
```

<br /> 

#### Launching job

```
# Launching a job

python
import requests
headers = {'Authorization': 'jacek.marzec@unimelb.edu.au abcd1234'}
payload = {'cancer_type': 'PAAD', 'title': 'example'}
r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
headers=headers,
files={'mutations': open('/data/example.txt', 'rb')}
,data=payload)
r.json()


# Getting a list with the identifiers of the jobs

headers = {'Authorization': 'jacek.marzec@unimelb.edu.au abcd1234'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
r.json()


# Accessing basic information about the job

headers = {'Authorization': 'jacek.marzec@unimelb.edu.au abcd1234'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/123456789', headers=headers)
r.json()
```

###### Note

The possible parameters for a job are:

* cancer_type - cancer type under analysis. This parameter is required and must be one of the values between brackets that you can find in the [analysis page](https://www.cancergenomeinterpreter.org/analysis).
* title - title for this run (optional)
* cnas - file containing the CNA data (more information about the formats in provided in the [FAQs page](https://www.cancergenomeinterpreter.org/faq#q13))
* mutations - file containing the mutation data (more information about the formats can be found in the [FAQs page](https://www.cancergenomeinterpreter.org/faq#q13))
* translocations - file containing the translocation data (more information about the formats in the [FAQs page](https://www.cancergenomeinterpreter.org/faq#q13))

The files (cnas, mutations and traslocations) are optional, but at least one is required.

<br /> 

#### Downloading results

Once the `job status` is indicated as `done` then once can download the results.

```
# Downloading results

headers = {'Authorization': 'jacek.marzec@unimelb.edu.au abcd1234'}
payload={'action':'download'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/123456789', headers=headers, params=payload)
with open('/data/example_cgi_analysis.zip', 'wb') as fd:
    fd.write(r._content)


# Uncompressing folder with results

unzip /data/example_cgi_analysis.zip -d /data/example_cgi_analysis
rm /data/example_cgi_analysis.zip


##### Downloading logs

headers = {'Authorization': 'jacek.marzec@unimelb.edu.au abcd1234'}
payload={'action':'logs'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/123456789', headers=headers, params=payload)
with open('/data/example_cgi_analysis/logs.txt', 'wb') as fd:
    fd.write(r._content)
```

The job and results can be deleted once there is no need to store the results on [CGI](https://www.cancergenomeinterpreter.org/api/v1) server.

```
##### Deleting a job
import requests
headers = {'Authorization': 'jacek.marzec@unimelb.edu.au abcd1234'}
r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/123456789', headers=headers)
r.json()
```

<br />

### Output

[CGI](https://www.cancergenomeinterpreter.org/api/v1) results folder contains four output files:

* `input01.tsv` - the converted MAF file that is used as an input for the analysis
* `mutation_analysis.tsv` - table containing annotations of validated oncogenic alterations and cancer driver genes
* `drug_prescription.tsv` - table with annotated alterations described as biomarkers for the selected tumor type
* `drug_prescription_bioactivities.tsv` - table with annotated gene-compound interactions

<br />

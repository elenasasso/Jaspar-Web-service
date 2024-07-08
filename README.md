# Scientific programming R project: Gene Classes Package


## Background
This project was developed as part of the Scientific programming course within the Master's degree in Bioinformatics for Computational Genomics.

## Authors
This project was developed by [Elena Sasso](https://github.com/elenasasso) .


## Project Overview
This project aims to create a RESTful Web Service that handles information about all the TF (Transcription Factor) binding motifs extracted from the 
JASPAR database (https://jaspar.elixir.no/). 
Each motif is defined by:
- a JASPAR motif ID
- the corresponding TF name
- the Position Frequency Matrix (PFM) of the motif.
- 
The web service provides the following functions to users:
- Obtain data for all stored motifs
- Obtain data for a single, specific stored motif
- Add a new motif to the database
- Update an existing motif in the database
- Delete an existing motif from the database
- Obtain the visual representation for a single, specific stored motif
- Get the sequence match score for each motif that has the same length of a DNA sequence submitted by the user
- Get the sequence match score for each motif that has the same length of an RNA sequence submitted by the user, after converting it into DNA
- Get the reverse translation of a protein sequence submitted by the user into all possible DNA sequences, to choose one of them to request the score for.

## How to Use This Repository
To start the web service, the internet connection is required to download the Jaspar database. Then it is sufficient to write from the command line the name of 
the script followed by the port number on which one desires to run it (or the -h to get the helper). First download the [JASPAR_WEB_SERVICE.py](JASPAR_WEB_SERVICE.py)
script and make sure you have python installed.

Example: 
(python JASPAR_WEB_SERVICE.py 5000)

(python JASPAR_WEB_SERVICE.py -h)

## Examples

To use the web service, after starting it on the command line with the above command, you can write in another window of the command line one of the following commands,
depending on what you like to get:

- To get all motifs: curl -i http://localhost:5000/Motifs/motif

- To get a single motif: curl -i http://localhost:5000/Motifs/motif/MA0004.1

- To add a new motif: curl -i -H "Content-type: application/json" -X POST -d "{\"motif_id\": \"MA1234.1\", \"TF_name\":\"name\", \"PFM\": {\"A\": [800, 807, 52, 61, 884, 851], \"C\": [68, 52, 29, 35, 22, 46], \"G\": [47, 44, 22, 28, 43, 42], \"T\": [85, 98, 898, 876, 51, 61]}}" http://localhost:5000/Motifs/motif

- To modify an existing motif: curl -i -H "Content-type: application/json" -X PUT -d "{\"TF_name\":\"name\", \"PFM\": {\"A\": [800, 807, 52, 61, 884, 851], \"C\": [68, 52, 29, 35, 22, 46], \"G\": [47, 44, 22, 28, 43, 42], \"T\": [85, 98, 898, 876, 51, 61]}}" http://localhost:5000/Motifs/motif/MA1234.1

- To delete an existing motif: curl -i -X DELETE http://localhost:5000/Motifs/motif/MA1234.1

- To get the scores for a given sequence: curl -i  http://localhost:5000/Motifs/ATGC

- To get the visual representation for a single motif: curl -i http://localhost:5000/Motifs/motif/graph/MA0004.1


## Contact
For any additional questions or feedback, please contact [Elena Sasso](mailto:elenasasso01@gmail.com).

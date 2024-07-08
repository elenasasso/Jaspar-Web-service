# Import all packages

import sys
import requests
from flask import Flask, jsonify, request
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os


# Here below there is the helper returned when the user write the -h parameter

helper = ("\n"
          "Nome script: JASPAR MOTIF MATCHER WEB SERVICE\n"
          "\n"
          "Positional arguments:\n"
          "\tPort number\tArguments to specify on which port run the service. This must be a number between 1024 and 65535\n"
          "\n"
          "Options:\n"
          "\t-h, --help\tShow this help message and exit"
          "\n")

# These lines of code below are to make sure that the user specify a right port to run the service or the -h parameter

try:
    # sys.argv are the arguments specified by the user on the command line, if they are less than 2 it means the user
    # only write the name of the script, but not the port number or the -h
    if len(sys.argv) < 2:
        raise ValueError("ERROR: It's necessary to specify a parameter, either the port number to run the web service or the -h for see the helper")
    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print(helper)
        sys.exit(0)
    elif sys.argv[1].isnumeric():
        port = int(sys.argv[1])
        # the port before 1023 are usually reserved
        if port < 1023 or port > 65535:
            raise ValueError("ERROR: port number doesn't exist or is not available")
    else:
        raise ValueError('ERROR: Second parameter not valid, specifiy a port number or -h to see the helper')

except ValueError as ve:
    print(ve)
    sys.exit(1)


welcome_msg = ("\n"
          "WELCOME IN THE JASPAR WEB SERVICE! Here you can see motifs (also the graphical version), modify them, add or delete them.\n"
          "In addition you can put a sequence and you will get a score for all Jaspar motifs of the same length.\n"
          "\n"
          "Try to use all the options with curl! See the examples below:\n"
          "\n"
          "- To get all motifs: curl -i http://localhost:5000/Motifs/motif\n"
          "\n"
          "- To get just a single motif: curl -i http://localhost:5000/Motifs/motif/MA0004.1\n"
          "\n"
          '- To add a new motif: curl -i -H "Content-type: application/json" -X POST -d "{\\"motif_id\\": \\"MA1234.1\\", \\"TF_name\\":\\"name\\", \\"PFM\\": {\\"A\\": [800, 807, 52, 61, 884, 851], \\"C\\": [68, 52, 29, 35, 22, 46], \\"G\\": [47, 44, 22, 28, 43, 42], \\"T\\": [85, 98, 898, 876, 51, 61]}}" http://localhost:5000/Motifs/motif\n'
          "\n"
          '- To modify an existing motif: curl -i -H "Content-type: application/json" -X PUT -d "{\\"TF_name\\":\\"name\\", \\"PFM\\": {\\"A\\": [800, 807, 52, 61, 884, 851], \\"C\\": [68, 52, 29, 35, 22, 46], \\"G\\": [47, 44, 22, 28, 43, 42], \\"T\\": [85, 98, 898, 876, 51, 61]}}" http://localhost:5000/Motifs/motif/MA1234.1\n'
          "\n"
          "- To delete an existing motif: curl -i -X DELETE http://localhost:5000/Motifs/motif/MA1234.1\n"
          "\n"
          "- To get the scores for a given sequence: curl -i  http://localhost:5000/Motifs/ATGC\n"
          "\n"
          "- To get the visual representation for a motif: curl -i http://localhost:5000/Motifs/motif/graph/MA0004.1\n"
          "\n"
          "Just be carefull to typos: don't forget any brackets or escape characters or you will get errors.\n"
          "The suggestion is to copy the commands in the examples above and modify just the part that you need without changing the format.\n"
          "Also make sure to have the internet connection."
          "\n"
          "Now get ready to start and have fun!"
          "\n")

print(welcome_msg)

# The following block of code (about 50 lines) is to create a local copy of the Jaspar database taking it
# directly from the web, without need to download it first, so there is no problem of path
# you just need the internet connection
# The final database will be a list of motifs, which will be dictionaries with 3 keys: motif_id, TF_name and PFM
# which also will be a dictionaries with 4 keys (the 4 nucleotides) and as values the frequencies

# URL of the Jaspar txt file
url = 'https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt'

response = requests.get(url)

if response.status_code == 200:
    with open('downloaded_file.txt', 'wb') as f:
        f.write(response.content)
else:
    print('Failed to download the file')

motifs_data = []

with open('downloaded_file.txt', 'r') as file:
    lines = file.readlines()

for line in lines:
    # if the line starts with '>' it means it's a new motif, it' a line with the motif_id and the TF_name
    # So I extract them and put them in a local variable (current_motif) to save data of the current motif
    # I also prepare the pfm as an empty dictionary
    # Then put all the data in the list of motifs
    if line.startswith('>'):
        current_motif = {}
        motif_id, tf_name = line.strip().split('\t')[0].split(' ')[0][1:], line.strip().split('\t')[1]
        current_motif['motif_id'] = motif_id
        current_motif['TF_name'] = tf_name
        current_motif['PFM'] = {'A': [], 'C': [], 'G': [], 'T': []}
        motifs_data.append(current_motif)
    else:
        # if it doesn't start with that symbol it means that it's the pfm, so I complete the values with the frequencies
        nucleotide, freqs_str = line.strip().split('[')[0].strip(), line.strip().split('[')[1].split(']')[0]
        current_motif['PFM'][nucleotide] = list(map(int, freqs_str.split()))


# The function below is to display a nicer output to the user when he asks for motifs (for example in the get method),
# instead of the dictionary
def format_all_motifs(motifs_data):
    formatted_output = ""
    for motif in motifs_data:
        formatted_output += f"Motif ID: {motif['motif_id']}\n"
        formatted_output += f"TF Name: {motif['TF_name']}\n"
        formatted_output += "PFM:\n"
        for nucleotide, freqs in motif['PFM'].items():
            formatted_output += f"{nucleotide}: {freqs}\n"
        formatted_output += "\n"
    return formatted_output


# The two function below are also to display nicer outputs, but they are needed with the jsonify objects that are
# returned to the user when curl something

def format_jsonify_output(json_output):
    formatted_output = ""
    for key, value in json_output.items():
        formatted_output += f"{key}\n{value}\n"
    return formatted_output


def jsonify_formatted(data, status=200):
    response = jsonify(data)
    response.status_code = status
    response.data = format_jsonify_output(data)
    return response


# Now I'm creating the web service: in the following mostly 150 lines of code I'm creating all the basic
# functionality to provide the user with a complete web service

app = Flask(__name__)


# Function to get all the motifs (calling the functions defined above to display them in a nicer way)
@app.route('/Motifs/motif', methods=['GET'])
def getAllMotifs():
    return jsonify_formatted({'Motifs:': format_all_motifs(motifs_data)})


# This function below is to check if a motif is really in the database
def check_motif_id_exists(motif_id):
    motif_ids = [motif['motif_id'] for motif in motifs_data]
    return motif_id in motif_ids


# Function to get just one specific motif
@app.route('/Motifs/motif/<motif_id>', methods=['GET'])
def getMotif(motif_id):
    try:
        # First I control if the motif is really in the database
        if not check_motif_id_exists(motif_id):
            raise ValueError('Motif ID not found')

        m = [motif for motif in motifs_data if (motif['motif_id'] == motif_id)]
        return jsonify_formatted({f'motif {motif_id}:': format_all_motifs(m)})

    except ValueError as ve:
        return jsonify_formatted({'ERROR:': str(ve)}), 400


# The following function is to check if the PFM has the right format
def validate_pfm(pfm):
    # Check if PFM data is a dictionary
    if not isinstance(pfm, dict):
        raise ValueError('PFM data is not a dictionary')

    # Check if keys are correct and in the correct order
    expected_keys = ['A', 'C', 'G', 'T']
    if list(pfm.keys()) != expected_keys:
        raise ValueError('Invalid PFM format: incorrect keys')

    # Check if values are lists of numbers
    length_of_list = None
    for key, value in pfm.items():
        if not isinstance(value, list) or not all(isinstance(num, (int, float)) for num in value):
            raise ValueError('Invalid PFM format: non-numeric values')

        if length_of_list is None:
            length_of_list = len(value)
        elif len(value) != length_of_list:
            raise ValueError('Invalid PFM format: list of nucleotides must be of the same dimension')


# Function to modify an already existing motif
@app.route('/Motifs/motif/<motif_id>', methods=['PUT'])
def updateMotif(motif_id):
    try:
        # First I control if the motif is really in the database
        if not check_motif_id_exists(motif_id):
            raise ValueError('Motif ID not found')

        m = [motif for motif in motifs_data if (motif['motif_id'] == motif_id)]
        data = request.json
        if 'TF_name' in data:
            tf_name = data.get('TF_name')
            m[0]['TF_name'] = tf_name

        if 'PFM' in data:
            pfm = data.get('PFM')

        # here I check if the pfm has the right format with the function defined above
            validate_pfm(pfm)

            m[0]['PFM'] = pfm

        # If everything okay return a statement of successful update
        return jsonify_formatted({f'updated motif {motif_id}': 'Check again the database to see the difference'})

    except ValueError as ve:
        return jsonify_formatted({'ERROR:': str(ve)}), 400


# function to post a new motif
@app.route('/Motifs/motif', methods=['POST'])
def CreateMotif():
    data = request.json
    motif_id = data.get('motif_id')
    tf_name = data.get('TF_name')
    pfm = data.get('PFM')

    try:
        # controls that all fields are present and not omitted in the post request
        if motif_id is None:
            raise ValueError('Motif ID is missing')
        if tf_name is None:
            raise ValueError('TF name is missing')
        if pfm is None:
            raise ValueError('PFM data is missing')

        # check if the motif already exists
        if check_motif_id_exists(motif_id):
            raise ValueError('Motif ID already exists')

        # control that the motif_id has the right pattern with the function provided in the re package
        # The right pattern must have MA as first two character, then 4 numbers.another number
        pattern = r'^MA\d{4}\.\d$'
        if not re.match(pattern, motif_id):
            raise ValueError('Invalid motif ID format')

        # here I check if the pfm has the right format with the function defined above
        validate_pfm(pfm)

        new_motif = {
            'motif_id': motif_id,
            'TF_name': tf_name,
            'PFM': pfm
        }
        motifs_data.append(new_motif)
        # If everything okay return a statement of successful post
        return jsonify_formatted({f'The motif {motif_id} was added': 'Check again the database to see it in the end!'})

    except ValueError as ve:
        return jsonify_formatted({'ERROR:': str(ve)}), 400

# Function to delete an already existing motif
@app.route('/Motifs/motif/<motif_id>', methods=['DELETE'])
def DeleteMotif(motif_id):
    try:
        # First I control if the motif is really in the database
        if not check_motif_id_exists(motif_id):
            raise ValueError('Motif ID not found')

        m = [motif for motif in motifs_data if (motif['motif_id'] == motif_id)]
        motifs_data.remove(m[0])
        return jsonify_formatted({'Deletion done': f'The motif {motif_id} was deleted'})

    except ValueError as ve:
        return jsonify_formatted({'ERROR:': str(ve)}), 400

# The following about 130 lines of codes are to create the additional function to get the score
# for all motifs of the same length, given a sequence
# I'm also giving the possibility to give an RNA sequence translating the U to T, and also if you try to give
# a protein sequence I'm going to give you all the possible DNA sequences for that protein sequence, so you can retry
# to get the score for one of the DNA sequences

# First function is just to display a nicer output for the scores
def format_frequenzeMotivi_output(frequenze_motivi):
    formatted_output = ""
    for key, value in frequenze_motivi.items():
        formatted_output += f"{key}\t{value}\n"
    return formatted_output


# This function is to normalize the numbers in the pfm so the scores are normalized
def normalizza_pfm(pfm):
    pfm_normalizzata = {}
    numero_colonne = len(list(pfm.values())[0])

    # Get the sum for each column
    somma_colonne = [sum(pfm[nucleotide][posizione] for nucleotide in pfm) for posizione in range(numero_colonne)]

    # Normalize numbers in the column
    for nucleotide in pfm:
        pfm_normalizzata[nucleotide] = [pfm[nucleotide][posizione] / somma_colonne[posizione] for posizione in
                                        range(numero_colonne)]

    return pfm_normalizzata


# Function to compute the score: I'm multiplying the frequencies for each position in the motif corresponding to the same letter in the sequence
def calcola_score(sequenza, pfm_normalizzata):
    score = 1
    for posizione, nucleotide in enumerate(sequenza):
        score *= pfm_normalizzata[nucleotide][posizione]
    return score


# function to translate the proteic sequence into all the DNA sequences
def dna_sequences_from_protein(protein_sequence, gencode):
    dna_sequences = []

    for amino_acid in protein_sequence:
        if amino_acid in gencode:
            dna_sequences.append(gencode[amino_acid])

    # generate all the possible sequences
    from itertools import product
    possible_dna_sequences = [''.join(seq) for seq in product(*dna_sequences)]

    return possible_dna_sequences

# genetic code
gencode = {'I': ['ATA', 'ATC', 'ATT'], 'M': ['ATG'], 'T': ['ACA', 'ACC', 'ACG', 'ACT'], 'N': ['AAC', 'AAT'], 'K': ['AAA', 'AAG'], 'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], 'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'P': ['CCA', 'CCC', 'CCG', 'CCT'], 'H': ['CAC', 'CAT'], 'Q': ['CAA', 'CAG'], 'V': ['GTA', 'GTC', 'GTG', 'GTT'], 'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'D': ['GAC', 'GAT'], 'E': ['GAA', 'GAG'], 'G': ['GGA', 'GGC', 'GGG', 'GGT'], 'F': ['TTC', 'TTT'], 'Y': ['TAC', 'TAT'], 'C': ['TGC', 'TGT'], 'W': ['TGG']}


# This is the function to get the score
@app.route('/Motifs/<sequenza>', methods=['GET'])
def getScore(sequenza):
    sequenza_originale = sequenza
    warning = None

    # control that in the sequence there are only letters
    try:
        if not sequenza.isalpha():
            raise ValueError('The sequence must include only letters')

        # Transforming all letters in uppercase
        sequenza = sequenza.upper()

        if any(letter in 'BJOZ' for letter in sequenza):
            raise ValueError('The sequence contains non-valid letters')

        if 'T' in sequenza and 'U' in sequenza:
            raise ValueError('The sequence contains both T and U. Please choose between DNA or RNA')

        # Now I'm distinguishing between different cases: DNA, RNA or protein
        # DNA
        if set(sequenza) <= set('ACGT'):
            dna_sequences = [sequenza]

        # RNA
        elif set(sequenza) <= set('ACGU'):
            warning = 'WARNING: The sequence contains U, it could be an RNA. I will replace U with T\n'
            sequenza = sequenza.replace('U', 'T')
            dna_sequences = [sequenza]

        # protein
        else:
            warning = 'WARNING: The sequence contains lots of letters, it could be a protein sequence. I will translate it back to DNA.\n'
            dna_sequences = dna_sequences_from_protein(sequenza, gencode)
            return jsonify_formatted({
                                f'{warning}Choose one options in the following DNA sequences list and then retry to get the score:': dna_sequences})

        lunghezza_sequenza = len(sequenza)

        # Find the maximum and minimum possible length of all motifs
        lunghezza_minima = min(len(motif['PFM']['A']) for motif in motifs_data)
        lunghezza_massima = max(len(motif['PFM']['A']) for motif in motifs_data)

        # Check if the given sequence is too long or too short
        if lunghezza_sequenza < lunghezza_minima:
            return jsonify_formatted({'ERROR': f"The sequence is too short, the minimum lenght of a motif is {lunghezza_minima}"})
        elif lunghezza_sequenza > lunghezza_massima:
            return jsonify_formatted({'ERROR': f"The sequence is too long, the maximum lenght of a motif is {lunghezza_massima}"})

        frequenze_motivi = {}

        # for each motif(after normalizing it) computing the score with the functions above
        for elemento in motifs_data:
            pfm = elemento['PFM']
            lunghezza_motivo = len(pfm['A'])  # all values have the same length
            for seq in dna_sequences:
                if lunghezza_motivo == len(seq):
                    pfm_normalizzata = normalizza_pfm(pfm)
                    score_motivo = calcola_score(seq, pfm_normalizzata)
                    frequenze_motivi[elemento['motif_id']] = score_motivo

        if not frequenze_motivi:
            raise ValueError('There are no motifs of the same length of your sequence')

        frequenze_motivi_ordinati = dict(sorted(frequenze_motivi.items(), key=lambda x: x[1], reverse=True))
        frequenze_motivi_bello = format_frequenzeMotivi_output(frequenze_motivi_ordinati)
        numero_motivi = len(frequenze_motivi_ordinati)

        if warning is not None:
            return jsonify_formatted({
                    f'{warning}Given sequence: {sequenza_originale}. Number of motifs analized: {numero_motivi}': f'Scores of motifs:\n{frequenze_motivi_bello}'})
        return jsonify_formatted({
                f'Given sequence: {sequenza_originale}. Number of motifs analized: {numero_motivi}': f'Score of motifs:\n{frequenze_motivi_bello}'})

    except ValueError as ve:
        return jsonify_formatted({'ERROR:': str(ve)}), 400


# Finally in the following lasts lines I'm doing an extra function to get a nice graph of the motif
# The graph will be save into an image in the same directory where the script is running
# and also the link to display it will be available after the curl


# Function to create the graph
def create_motif_graph(motif):
    pfm = motif['PFM']
    pfm_norm = normalizza_pfm(pfm)
    positions = np.arange(len(pfm_norm['A']))  # Posizioni nel motivo

    fig, ax = plt.subplots()
    width = 0.35

    # creation of a bar for each nucelotide
    bars = []
    bottom = np.zeros(len(pfm_norm['A']))  # initial position of the bars
    nucleotides = ['A', 'C', 'G', 'T']
    for nucleotide in nucleotides:
        bars.append(ax.bar(positions, pfm_norm[nucleotide], width, bottom=bottom, label=nucleotide))
        bottom += np.array(pfm_norm[nucleotide])  # update the position for the next bar

    # Add of labels and legend
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Position')
    ax.set_title('Motif Frequencies')
    ax.set_xticks(positions)
    ax.set_xticklabels(positions + 1)
    ax.legend()

    ax.text(0.5, 1.08, f'Motif ID: {motif["motif_id"]}', horizontalalignment='center', transform=ax.transAxes)

    # Save the graph as an image
    filename = f'motif_graph_{motif["motif_id"]}.png'
    plt.savefig(filename)
    plt.close()

    return filename


def get_image_link(filename):
    abs_path = os.path.abspath(filename)  # get the path of the image
    return f'file://{abs_path}'


# Function to obtain the graph of a single motif
@app.route('/Motifs/motif/graph/<motif_id>', methods=['GET'])
def getMotifGraph(motif_id):
    try:
        # First I control if the motif is really in the database
        if not check_motif_id_exists(motif_id):
            raise ValueError('Motif ID not found')

        motif = [motif for motif in motifs_data if motif['motif_id'] == motif_id][0]

        # Create the graph for the given motif
        motif_graph_filename = create_motif_graph(motif)
        image_link = get_image_link(motif_graph_filename)

        # return a message with the link of the image
        return jsonify({'image_link': image_link}), 200

    except ValueError as ve:
        return jsonify_formatted({'ERROR:': str(ve)}), 400


# This is just to catch errors deriving from some typos in the request of the users like forgotting a curly
# bracket at the end of the pfm during a post or a put request
@app.errorhandler(400)
def handle_bad_request(error):
    return jsonify_formatted({'ERROR': 'Bad request. Please check your input.'}), 400


@app.errorhandler(404)
def handle_not_found(error):
    return jsonify_formatted({'ERROR': 'Not found. Please check your input.'}), 404


if __name__ == '__main__':
    app.run(port=port)


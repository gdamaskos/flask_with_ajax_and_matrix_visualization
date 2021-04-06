import asyncio
import logging
import requests
import sys
import traceback
import urllib.request, urllib.error, urllib.parse

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from flask import (json, render_template,
                   request, send_from_directory,
                   send_file, url_for, redirect, jsonify)

from ccd import (app, sequence, dna_finder_async)
from ccd.custom_exceptions import (NoNeedToAlignError, NoIsoformsPresentError,
                                   DNASearchTimeoutError,
                                   IdNotFoundError, IGiveUpError)
from ccd.services import predict


_log = logging.getLogger(__name__)
ccd_sequences = SeqRecord(Seq(''))



@app.route('/')
def index():
    return render_template('index.html')


@app.route('/tutorial')
def tutorial():
    return render_template('tutorial.html')


@app.route('/credits')
def credits():
    return render_template('credits.html')


@app.route('/help')
def help():
    return render_template('help.html')


@app.route('/contact')
def contact():
    return render_template('contact.html')


@app.route('/search_uniprot/<uniprot_accession_or_entry_name>', methods=['GET'])
def search_uniprot(uniprot_accession_or_entry_name):
    _log.debug('Searching Uniprot..')
    try:
        loop = asyncio.new_event_loop() #get_event_loop does not seem to work here
        asyncio.set_event_loop(loop)
        search_object = dna_finder_async.DatabaseCrawler(uniprot_accession_or_entry_name)
        try:
            protein = search_object.search() #if this does not complete in timeout seconds, handler() will raise an Exception
        except DNASearchTimeoutError:
            _log.error('Downloading DNA entries from NCBI took too long')
            return jsonify({'status': 'ERROR', 'reason': 'Failed to fetch crossreferences from DNA databases.'})
        
        data = protein.serialize_protein()
        try: #checking if we can align isoforms
            aligned_isoforms = [s.serialize() for s in protein.align_isoforms()]
            data.update({'aligned_isoforms': aligned_isoforms})
        except (NoNeedToAlignError, NoIsoformsPresentError):
            pass
        return jsonify(data)
    except IdNotFoundError:
        _log.error('ID not found.')
        return jsonify({'status': 'ERROR', 'reason': 'ID not found.'})
    except IGiveUpError:
        _log.error('Unable to map entry to DNA.')
        return jsonify({'status': 'ERROR', 'reason': 'Unable to map entry to DNA.'})
    except urllib.error.HTTPError as err:
        if err.code in range(500, 506):
            msg = 'ID "{}" does not seem to be a valid uniprot accession number or identifier.'.format(uniprot_accession_or_entry_name)
            _log.error(err)
            return jsonify({'status': 'ERROR', 'reason': msg})
        elif err.code in range(400, 451):
            msg = 'Problem communicating with the remote databases, please try again later.'
            _log.error(err)
            return jsonify({'status': 'ERROR', 'reason': msg})
    except requests.HTTPError as e:
        if e.response.status_code in range(500, 506):
            msg = 'ID "{}" does not seem to be a valid uniprot accession number or identifier.'.format(uniprot_accession_or_entry_name)
            _log.error(e)
            return jsonify({'status': 'ERROR', 'reason': msg})
        elif e.response.status_code in range(400, 451):
            msg = 'Problem communicating with the remote databases, please try again later.'
            _log.error(e)
            return jsonify({'status': 'ERROR', 'reason': msg})
    except Exception as e:
        _log.error(e)
        traceback.print_exc(file=sys.stdout) 
        return jsonify({'status': 'ERROR', 'reason': 'Unknown error.'})


@app.route('/translate', methods=['POST'])
def translate():
    ccd_sequences.letter_annotations = {}
    dna_seq = request.values.get('DNAseq')
    if sequence.dna_validator(dna_seq):
        # If given, remove stop codon
        if dna_seq[-3:] in ['TAG', 'TGA', 'TAA']:
            dna_seq = dna_seq[:-3]
        # Amino acids sequence is needed for all predictions
        aa_seq_obj = sequence.translate(dna_seq)
        # Add to sequence record object the amino acid sequence and
        # the open reading frame DNA (in triplets)
        _log.debug('Adding data to sequence record biopython object..')
        ccd_sequences.seq = aa_seq_obj
        orf_triplets = []
        for start in range(0, len(dna_seq), 3):
            orf_triplets.append(dna_seq[start: start + 3])
        ccd_sequences.letter_annotations["orf_triplets"] = orf_triplets
        return json.dumps({'status': 'OK', 'AAseq': str(ccd_sequences.seq)})
    else:
        return json.dumps({'status': 'VALIDATE_ERROR'})


@app.route('/prediction/<service_name>', methods=['GET'])
def prediction(service_name):
    _log.debug('Getting predictions..')
    seq = str(ccd_sequences.seq)
    try:
        prediction = predict.service(service_name, seq)
        ccd_sequences.letter_annotations[prediction[0]] = prediction[1]
        return json.dumps({'status': 'OK',
                           'prediction': prediction,
                           'error': 'NO'})
    except TypeError as t:
        _log.debug('Something went wrong when creating Seqrecord annotations')
        _log.debug(f'prediction is of length {len(prediction[1])} and is {prediction[1]}')
        _log.debug(f'sequence is of length {len(seq)}')
        return 'Prediction Error'
    except Exception as e:
        traceback.print_exc(file=sys.stdout) 
        error_string = '{} {}'.format(service_name, str(e))
        _log.error(error_string)
        return json.dumps({'status': 'OK',
                           'prediction': service_name,
                           'error': 'YES'})
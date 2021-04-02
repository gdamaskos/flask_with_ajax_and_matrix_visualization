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

from ccd import (app, sequence, vector_sequence,
                 orf_finder, dna_finder_async,
                 homolog_finder, celery_tasks)
from ccd.custom_exceptions import (NoNeedToAlignError, NoIsoformsPresentError,
                                   DNASearchTimeoutError,
                                   IdNotFoundError, IGiveUpError)
from ccd.services import predict




_log = logging.getLogger(__name__)
ccd_sequences = SeqRecord(Seq(''))
evo = homolog_finder.Evolutionator()
evo.searched = False
timeout = 45


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


@app.route('/get_vectors', methods=['GET'])
def get_vectors():
    try:
        return send_from_directory('static', 'vectors_enzymes.json')
    except Exception as e:
        traceback.print_exc(file=sys.stdout) 
        _log.error(str(e))
        return str(e)


@app.route('/evolutionator', methods=['POST'])
def evolutionator():
    print(f'evo.searched is now {evo.searched}')
    uniprot_id = request.values.get('uniprot_id')
    aa_seq = request.values.get('AAseq')
    #uniprot_id = request.values.get('uniprot_id')
    flags_dictionary = json.loads(request.values.get('flags_dictionary'))
    try:
        if uniprot_id and not evo.searched:
            evo.fill_data(uniprot_id, aa_seq, '')
            evo.searched = True
            print(f'evo.searched is now {evo.searched}')
        else:
            evo.fill_data('', aa_seq, '')
        sequences = evo.search()
        flags, serialized_alignment, hits_by_flag = evo.get_default_alignment(
                                                    flags_dictionary, sequences)
        return json.dumps({'status': 'OK',
                           'multiple_alignment': serialized_alignment,
                           'alignment_flags': flags,
                           'hits_by_flag': hits_by_flag})
    except IdNotFoundError as e:
        msg = 'Could not complete the multiple alignment. ID not found.'
        _log.error('{} {} {}'.format(msg, str(e), type(e)))
        _log.error(msg)
        return json.dumps({'status': 'ERROR', 'reason': msg})
    except IGiveUpError:
        msg = 'Could not complete the multiple alignment. Unable to map entry to DNA.'
        _log.error(msg)
        return json.dumps({'status': 'ERROR', 'reason': msg})
    except urllib.error.HTTPError as err:
        if err.code in range(500, 506):
            msg = 'Could not complete the multiple alignment. Remote server error.'
            _log.error('{} {}'.format(msg, err))
            return json.dumps({'status': 'ERROR', 'reason': msg})
        elif err.code in range(400, 451):
            msg = 'Could not complete the multiple alignment. Problem communicating with the remote databases, please try again later.'
            _log.error('{} {}'.format(msg, err))
            return json.dumps({'status': 'ERROR', 'reason': msg})
    except requests.HTTPError as e:
        if e.response.status_code in range(500, 506):
            msg = 'Could not complete the multiple alignment. Remote server error.'
            _log.error('{} {}'.format(msg, e))
            return json.dumps({'status': 'ERROR', 'reason': msg})
        elif e.response.status_code in range(400, 451):
            msg = 'Could not complete the multiple alignment. Problem communicating with the remote databases, please try again later.'
            _log.error('{} {}'.format(msg, e))
            return json.dumps({'status': 'ERROR', 'reason': msg})
    except Exception as e:
        error_string = 'Could not complete the multiple alignment. Unknown error.'
        traceback.print_exc(file=sys.stdout)
        _log.error('{} {} {}'.format(error_string, str(e), type(e)))
        return json.dumps({'status': 'ERROR', 'reason': error_string})
    
@app.route('/realign', methods=['POST'])
def realign():
    flags_dictionary = json.loads(request.values.get('flags_dictionary'))
    try:
        flags, serialized_alignment = evo.realign(flags_dictionary, evo.cached_sequences)
        return json.dumps({'status': 'OK', 
                           'realignment': serialized_alignment,
                           'realignment_flags': flags})
    except IdNotFoundError:
        msg = 'Could not realign. ID not found.'
        _log.error(msg)
        return json.dumps({'status': 'ERROR', 'reason': msg})
    except IGiveUpError:
        msg = 'Could not realign. Unable to map entry to DNA.'
        _log.error(msg)
        return json.dumps({'status': 'ERROR', 'reason': msg})
    except urllib.error.HTTPError as err:
        if err.code in range(500, 506):
            msg = 'Could not realign. Remote server error.'
            _log.error('{} {}'.format(msg, err))
            return json.dumps({'status': 'ERROR', 'reason': msg})
        elif err.code in range(400, 451):
            msg = 'Could not realign. Problem communicating with the remote databases, please try again later.'
            _log.error('{} {}'.format(msg, err))
            return json.dumps({'status': 'ERROR', 'reason': msg})
    except requests.HTTPError as e:
        if e.response.status_code in range(500, 506):
            msg = 'Could not realign. Remote server error.'
            _log.error('{} {}'.format(msg, e))
            return json.dumps({'status': 'ERROR', 'reason': msg})
        elif e.response.status_code in range(400, 451):
            msg = 'Could not realign. Problem communicating with the remote databases, please try again later.'
            _log.error('{} {}'.format(msg, e))
            return json.dumps({'status': 'ERROR', 'reason': msg})
    except Exception as e:
        error_string = 'Could not realign.'
        traceback.print_exc(file=sys.stdout) 
        _log.error('{} {}'.format(error_string, str(e)))
        return json.dumps({'status': 'ERROR', 'reason': error_string})


@app.route('/malignment_fasta', methods = ['GET'])
def malignment_fasta():
    flags_dictionary = json.loads(request.values.get('flags_dictionary'))
    new_alignment = evo.realign(flags_dictionary, evo.cached_sequences)
    strIO = evo.get_fasta(new_alignment)
    return send_file(strIO,
                     attachment_filename='ccd_multiple_alignment.fasta',
                     as_attachment=True)

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
#     except DNASearchTimeoutError:
#         _log.error('Downloading DNA entries from NCBI took too long')
#         return jsonify({'status': 'ERROR', 'reason': 'Failed to fetch crossreferences from DNA databases.'})
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


@app.route('/primersnpeptides', methods=['POST'])
def primersnpeptides():
    starts = request.values.get('starts')
    starts = json.loads(starts)

    stops = request.values.get('stops')
    stops = json.loads(stops)

    method_value = request.values.get('method_value')
    method_value = json.loads(method_value)

    forward_overhang = request.values.get('FWoverhang')
    reverse_overhang = request.values.get('RVoverhang')
    if sequence.over_validator(forward_overhang) and \
            sequence.over_validator(reverse_overhang):
        _log.debug('Building primers and validating starts/stops..')
        dna_seq = ''.join(ccd_sequences.letter_annotations["orf_triplets"])
        if method_value[0] == 0:  # use melting temperature
            primers, errors, vstarts, vstops = \
                    sequence.build_primers_tm(dna_seq, starts, stops,
                                              forward_overhang,
                                              reverse_overhang,
                                              int(method_value[1]))
            MWs, pIs, epsilons = sequence.protparams(str(ccd_sequences.seq),
                                                     vstarts, vstops)
        if method_value[0] == 1:  # use number of bases
            primers, errors, vstarts, vstops = \
                    sequence.build_primers_bs(dna_seq, starts, stops,
                                              forward_overhang,
                                              reverse_overhang,
                                              int(method_value[1]))
            MWs, pIs, epsilons = sequence.protparams(str(ccd_sequences.seq),
                                                     vstarts, vstops)
        return json.dumps({'status': 'OK',
                           'primers': primers,
                           'vstarts': vstarts,
                           'vstops': vstops,
                           'MWs': MWs,
                           'pIs': pIs,
                           'epsilons': epsilons,
                           'errors': errors})
    else:
        return json.dumps({'status': 'VALIDATE_ERROR'})


@app.route('/get_separated_vectors', methods = ['GET'])
def get_separated_vectors():
    try:
        vnames = vector_sequence.vector_names()
        return json.dumps({'status': 'OK', 'vnames': vnames})
    except Exception as e:
        traceback.print_exc(file=sys.stdout) 
        _log.error(str(e))
        return json.dumps({'status': 'VECTORS_ERROR'})


@app.route('/mapmake', methods = ['GET'])
def mapmake():
    vname = request.values.get('vname')
    pname = request.values.get('pname')
    pmfname = request.values.get('pmfname')
    pstarts = request.args.getlist('pstarts[]')
    pstops = request.args.getlist('pstops[]')
    
    dna_seq = ''.join(ccd_sequences.letter_annotations["orf_triplets"])
    inserts_dna_seq = vector_sequence.poly2dna_seq(dna_seq, pname, pstarts, pstops)
    zbuff = vector_sequence.make_plasmid_maps(inserts_dna_seq, vname)

    return send_file(zbuff,
                     attachment_filename=pmfname,
                     as_attachment=True)


@app.route('/crysol_prediction/<service_name>', methods=['POST'])
def crysol_prediction(service_name):
    pp_seqs = request.values.get('selected_polypeptides_sequences')
    pp_seqs = json.loads(pp_seqs)

    _log.debug('Getting crysol predictions..')
    crysol_task = celery_tasks.crysol_service.delay(service_name, pp_seqs)
    return jsonify({}), 202, {'Location': url_for('taskstatus', _external=True, _scheme='https',
                                                  task_id=crysol_task.id)}
    
    
@app.route('/status/<task_id>', methods=['GET'])
def taskstatus(task_id):
    task = celery_tasks.crysol_service.AsyncResult(task_id)
    response = {}
    
    if task.state == 'SUCCESS':
        response = {
            'state': task.state,
            'result': task.get(propagate=False),
            'info': str(task.info)
        }
    # defensive programming below, because success state was passing the first check sometimes
    elif task.state == 'FAILURE' or task.state == 'PENDING' or \
            task.state == 'REVOKED' or task.state == 'RETRY' or task.state == 'RECEIVED':
        # other state or something went wrong in the background job
        response = {
            'state': task.state,
            'geom': 'nullis'
        }
    return jsonify(response)
 
 
@app.route('/revokeall', methods=['GET'])
def alltasksrevoke():
    task_ids_list = request.args.getlist('task_ids_list[]')
    for task_id in task_ids_list:
        task = celery_tasks.crysol_service.AsyncResult(task_id)
        task.revoke(terminate=True)
    return jsonify({}), 202, {'status': 'OK'}  
     
    
@app.route('/tag_cleave', methods=['POST'])
def tag_cleave():
    polypeptides = request.values.get('polypeptides')
    polypeptides = json.loads(polypeptides)
    pstarts = [row[2] for row in polypeptides]
    pstops = [row[3] for row in polypeptides]
    vname = request.values.get('vname')
    pname = request.values.get('pname')
    
    dna_seq = ''.join(ccd_sequences.letter_annotations["orf_triplets"])
    inserts_dna_seq = vector_sequence.poly2dna_seq(dna_seq, pname, pstarts, pstops)
    tagged_cleaved = orf_finder.tag_cleave_pps(inserts_dna_seq, vname)   
    
    return json.dumps({'tagged_cleaved': tagged_cleaved})


# A bunch of redirects from the old ccd
# Obtained from 404s in the logs
@app.route('/welcome.html')
def redirect_welcome():
    return redirect(url_for('index'), 301)


@app.route('/Welcome.html')
def redirect_welcome_capital():
    return redirect(url_for('index'), 301)


@app.route('/design.html')
def redirect_design():
    return redirect(url_for('index'), 301)


@app.route('/Design.html')
def redirect_design_capital():
    return redirect(url_for('index'), 301)


@app.route('/help.html')
def redirect_help():
    return redirect(url_for('help'), 301)


@app.route('/Help.html')
def redirect_help_capital():
    return redirect(url_for('help'), 301)


@app.route('/start_ccd.html')
def redirect_start_ccd():
    return redirect(url_for('index'), 301)


@app.route('/PCR_oligos.html')
def redirect_PCR_oligos():
    return redirect(url_for('index'), 301)


@app.route('/Inspect.html')
def redirect_inspect():
    return redirect(url_for('index'), 301)


@app.route('/FAQ.html')
def redirect_faq():
    return redirect(url_for('help'), 301)


@app.route('/Input_DNA.html')
def redirect_input_DNA():
    return redirect(url_for('index'), 301)


@app.route('/Tutorial.html')
def redirect_tutorial():
    return redirect(url_for('tutorial'), 301)


@app.route('/Meta-actions.html')
def redirect_meta_actions():
    return redirect(url_for('index'), 301)

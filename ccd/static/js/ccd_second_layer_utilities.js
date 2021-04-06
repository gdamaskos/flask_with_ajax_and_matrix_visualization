// clear time!
'use strict';

function clear_uniprot_search() {
    $('#uniprot-id-search-alerts').html('');
    $('#uniprot-id-search-results').html('');
    $('.uniprot-search-progress-bar').css('display', 'none');
    $('#dna-orf').html('');

    if (CCD.unisearch_ajax) {
        CCD.unisearch_ajax.abort();
    }
}


function clear_all() {
    CCD.query_str = null;
    CCD.predictions = null;
    CCD.gapped_predictions = null;
    $('#predictions-alerts-placeholder').html('');
    $('#predictions-info-placeholder').html('');
    $('sequences-labels').html('');
    $('sequences-viewer').html('');
    $('aa-label').html('');
    $('aa-viewer').html('');
}



// hide time
function hide_byclass(classname) {
    $('.' + classname).each(function() {
        $(this).css('display', 'none');
    });
}


// show time!
function show_byclass(classname) {
    $('.' + classname).each(function() {
        $(this).css('display', 'block');
    });
}


function validate_search_uniprot_input(uniprot_id) {

    if (uniprot_id === null || uniprot_id === '') {
        $('#uniprot-id-search-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
            '<span>No ID to search for. Please provide an ID.</span></div>');
        return false;
    }

    // if the given uniprot id is not in accession format and not in entry name format, then notify user for the error.
    if (uniprot_id.search(/^[a-zA-Z0-9]{1,10}_[a-zA-Z0-9]{1,5}$/g) === -1 &&
        uniprot_id.search(/^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$/gi) === -1) {

        $('#uniprot-id-search-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
            '<span>The given ID is not a valid uniprot accession number or entry name.</span></div>');
        return false;
    }
}


function draw_aa() {
    var AAseq_span = '';
    var indexes_string = '';
    var curr_length = 0;
    var i;
    var bgcolor = '';
    var index_str = '';

    for (i = 0; i < CCD.AAseq.length; i += 1) {
        bgcolor = '';
        index_str = (i + 1).toString();
        // if index can be found for this amino acid number, color
        if (CCD.starts.indexOf(index_str) !== -1) {
            bgcolor = 'lawngreen';
        }
        if (CCD.stops.indexOf(index_str) !== -1) {
            bgcolor = 'tomato';
        }
        if (CCD.starts.indexOf(index_str) !== -1 && CCD.stops.indexOf(index_str) !== -1) {
            bgcolor = 'yellow';
        }
        AAseq_span += '<span class="amino-acid" style="cursor:pointer; background-color:' + bgcolor + ';" onclick="onclickAA(event)" data-toggle="tooltip" title="residue No ' +
            (i + 1) + '" name="' + (i + 1) + '">' + CCD.AAseq.substring(i, i + 1) + '</span>';
    }
    $('#aa-label').html('<span><font color="red">Query Seq</font></span><br><br>');
    $('#aa-viewer').html(AAseq_span + '<br>');
    // here we add the indexes as a string
    while (curr_length < CCD.AAseq.length - 1) {
        if (curr_length % 10 === 9) {
            indexes_string += '|' + (curr_length + 1).toString();
        } else {
            indexes_string += ' ';
        }
        curr_length = indexes_string.length;
    }
    $('#aa-viewer').append('<span style="white-space:pre;">' + indexes_string + '</span><br>');
}


function draw_predictions() {
    $('#sequences-labels').html('');
    $('#sequences-viewer').html('');
    //draw errors first
    var formatted_errors = '';
    var i;

    if (CCD.pred_errors.length !== 0) {
        for (i = 0; i < CCD.pred_errors.length; i += 1) {
            formatted_errors += CCD.pred_errors[i] + ' error<br>';
        }
        $('#predictions-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a><span>' +
            formatted_errors + '</span></div>');
    }
    $('.service-label').each(function() {
        var prediction_drawn = false;

        if ($(this).children('input').is(':checked') === true) {
            for (i = 0; i < CCD.pred_errors.length; i += 1) {
                if (CCD.pred_errors[i] === $(this).text().substr(3)) { // if there is an error in the name of this prediction name..
                    return true; // continues the .each function to the next element!
                }
            }
            if (CCD.predictions) {
                // recreate gapped_predictions in light of the new data in predictions, a new line
                var alignment_visibility = $('#alignment-toggler').is(':checked');
                if (CCD.query_str && alignment_visibility === true) { // null value is falsy, and means that alignment has not arrived
                    CCD.gapped_predictions = new Array(CCD.predictions.length);
                    for (i = 0; i < CCD.gapped_predictions.length; i += 1) {
                        CCD.gapped_predictions[i] = CCD.predictions[i].slice();
                    }
                    var gaps_counter = 0;
                    var gaps_start = 0;
                    var gaps_spot = [];
                    var position = 0;
                    var gaps_insert = '';
                    var new_prediction = '';
                    var new_indexes_string = '';
                    var indexes_string = '';
                    var curr_length = 0;
                    var AAseq_span = '';
                    var residue_counter = 1;
                    var j, k, n;
                    var bgcolor = '';
                    var index_str = '';

                    // performing: position finding and counting of the gaps
                    for (k = 0; k < CCD.query_str.length; k += 1) {
                        if (CCD.query_str[k] === ' ') {
                            if (gaps_counter === 0) {
                                gaps_start = k;
                            }
                            gaps_counter += 1;
                        } else {
                            // if counter is not zero, push index and counter 
                            if (gaps_counter !== 0) {
                                gaps_spot.push([gaps_start, gaps_counter]);
                                gaps_counter = 0;
                            }
                        }
                    }
                    if (gaps_counter !== 0) {
                        gaps_spot.push([gaps_start, gaps_counter]);
                    }

                    // draw again aa sequence
                    AAseq_span = '';
                    for (j = 0; j < CCD.query_str.length; j += 1) {
                        if (CCD.query_str[j] !== ' ') {
                            bgcolor = '';
                            index_str = residue_counter.toString();
                            // if index can be found for this amino acid number, color
                            if (CCD.starts.indexOf(index_str) !== -1) {
                                bgcolor = 'lawngreen';
                            }
                            if (CCD.stops.indexOf(index_str) !== -1) {
                                bgcolor = 'tomato';
                            }
                            if (CCD.starts.indexOf(index_str) !== -1 && CCD.stops.indexOf(index_str) !== -1) {
                                bgcolor = 'yellow';
                            }
                            AAseq_span += '<span class="amino-acid" style="cursor:pointer; background-color:' + bgcolor + ';" onclick="onclickAA(event)" data-toggle="tooltip" title="residue No ' +
                                residue_counter + '" name="' + residue_counter + '">' + CCD.query_str.substring(j, j + 1) + '</span>';
                            residue_counter += 1;
                        } else {
                            AAseq_span += '<span class="alignment-gap">&nbsp</span>';
                        }
                    }
                    $('#aa-label').html('<font color="red">Query Seq</font>' + '<br><br>');
                    $('#aa-viewer').html('<span style="white-space:pre;">' + AAseq_span + '</span><br>');
                    // here we add the indexes as a string
                    while (curr_length < CCD.AAseq.length - 1) {
                        if (curr_length % 10 === 9) {
                            indexes_string += '|' + (curr_length + 1).toString();
                        } else {
                            indexes_string += ' ';
                        }
                        curr_length = indexes_string.length;
                    }

                    new_indexes_string = indexes_string;
                    for (k = 0; k < gaps_spot.length; k += 1) {
                        position = gaps_spot[k][0];
                        gaps_insert = Array(gaps_spot[k][1] + 1).join(' ');
                        new_indexes_string = [new_indexes_string.slice(0, position), gaps_insert, new_indexes_string.slice(position)].join('');
                    }

                    $('#aa-viewer').append('<span style="white-space:pre;">' + new_indexes_string + '</span><br>');

                    // gapped_predictions
                    var gaps_insert = '';
                    //gaps_spot = Array(Array(pos1, length1), Array(pos2,length2), ...) where pos is the start of the gap and length is the length
                    //for each gapped position, we make a string of spaces of the correct length
                    for (k = 0; k < gaps_spot.length; k += 1) {
                        position = gaps_spot[k][0];
                        gaps_insert = Array(gaps_spot[k][1] + 1).join(' ');
                        total_gap_length += gaps_spot[k][1];
                        // and we add that string of spaces to each prediction
                        for (n = 0; n < CCD.gapped_predictions.length; n += 1) { // for every dotted_prediction element
                            CCD.gapped_predictions[n][1] = [CCD.gapped_predictions[n][1].slice(0, position), gaps_insert,
                                CCD.gapped_predictions[n][1].slice(position)
                            ].join('');
                            if (CCD.gapped_predictions[n][2]) { // some predictions have an extra element, which also needs to be gapped
                                CCD.gapped_predictions[n][2] = [CCD.gapped_predictions[n][2].slice(0, position),
                                    gaps_insert, CCD.gapped_predictions[n][2].slice(position)
                                ].join('');
                            }
                        }
                    }

                    //update the positions of start and stop points in PDB_95 predictions to account for gaps
                    //PDB_95 -> CCD.gapped_predictions[3] = JSON object {"pos1" = Array(structure1, structure2, ...), ...}
                    //since I cannot understand how gaps_spot is calculated by George and how those numbers map to the original
                    //sequence position, I will shamelessly exploit the fact that CCD.gapped_predictions[3] is ordered ascending
                    var total_gap_length;
                    var new_pdb_points = {};
                    var old_pos, new_pos;
                    for (var n = 0; n < CCD.gapped_predictions.length; n += 1) {
                        if (CCD.gapped_predictions[n][0] === 'PDB_95') {
                            //first we need to map gaps_spot back to the sequence's number. 
                            var seq = CCD.gapped_predictions[n][1];
                            var structures = JSON.parse(CCD.gapped_predictions[n][3]);
                            var positions = Object.keys(structures)
                            var counter = 0;
                            for (var i = 0; i <= seq.length; i += 1) {
                                if (seq[i] === '>' || seq[i] === '<' || seq[i] === String.fromCharCode('9674')) {
                                    new_pdb_points[i] = structures[positions[counter]];
                                    counter += 1;
                                }
                            }
                            CCD.gapped_predictions[n][3] = JSON.stringify(new_pdb_points);
                        }
                    }



                } else {
                    draw_aa();
                    CCD.gapped_predictions = CCD.predictions;
                }
            }
            // loop through the predictions that have arrived and gapped, if there draw it
            for (i = 0; i < CCD.gapped_predictions.length; i += 1) {
                if (CCD.gapped_predictions[i][0] === $(this).text().substr(3)) { // if prediction has arrived
                    switch ($(this).text().substr(3)) {
                        case 'PREDATOR':
                            draw_he(i);
                            break;
                        case 'IUPRED':
                        case 'GLOBPLOT':
                            draw_domains(i);
                            break;
                        case 'SMART':
                            draw_smart(i);
                            break;
                        case 'COILS':
                            draw_coils(i);
                            break;
                    }
                    prediction_drawn = true;
                }
            }
            // if nothing drawn, draw wait since it is checked
            if (prediction_drawn === false) {
                draw_wait($(this).text().substr(3));
            }
            // scroll bar check and pad
            var sbHeight = window.innerWidth - document.documentElement.clientWidth;
            if (window.innerWidth > document.documentElement.clientWidth) {
                document.getElementById("sequences-labels").style.paddingBottom = sbHeight + "px";
            } else {
                document.getElementById("sequences-labels").style.paddingBottom = "0px";
            }
        }
    });
}


function validateDNA(seq) {
    //Based on: http://www.blopig.com/blog/2013/03/a-javascript-function-to-validate-fasta-sequences/
    var i;
    var match;

    // clear previous alert messages
    $('#input-alerts').html('');
    // remove numbers
    seq = seq.replace(/\d/g, '');
    // remove trailing spaces
    seq = seq.trim();
    if (seq === null || seq === '') {
        $('#input-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
            '<span>No DNA to submit. Please input DNA sequence.</span></div>');
        return false;
    }
    // split on newlines
    var lines = seq.split('\n');
    // check for lines that begin with semicolon or grater symbol and remove them - remove fasta headers
    i = 0;
    while (i < lines.length) {
        if (lines[i][0] === ';' || lines[i][0] === '>') {
            lines.splice(i, 1); // remove one line, the current one 
            i = 0; // go back to the start
        } else {
            i += 1;
        }
    }
    // join the array back into a single string without newlines and
    // trailing or leading spaces
    seq = lines.join('').trim();
    // loose spaces in between and upper case
    seq = seq.replace(/[\s]/g, "").toUpperCase();
    // Check if something is still there
    if (seq === null || seq === '') {
        $('#input-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
            '<span>No DNA to submit after removing comments and headers. Please input DNA sequence.</span></div>');
        return false;
    }

    var warnings = '';
    if (seq.search(/^ATG/g) === -1) {
        // The seq string lacks the start codon
        warnings += 'Your ORF should start with an ATG.<br>';
    }
    if (seq.search(/(TGA|TAA|TAG)$/g) === -1) {
        // The seq string lacks stop codon
        warnings += 'Your ORF should end with a stop codon (TGA, TAA, TAG).<br>';
        if (seq.length < 93) {
            warnings += 'Your ORF should be at least 93 bp long.<br>';
        }
    } else if (seq.length < 90) {
        //The seq string contains should code at least 30 amino acids for some software to run smoothly, not taking into account stop codon that will be removed later from backend
        warnings += 'Your ORF should be at least 93 bp long.<br>';
    }

    // Search for charaters that are not G, A, T or C.
    if (seq.search(/^[ATGC]+$/g) === -1) {
        //The seq string contains non-DNA characters
        var nonDNA = /[^ATGC]/g;
        var matches = '';
        match = nonDNA.exec(seq);

        while (match) {
            matches += (match.index + 1).toString() + ', ';
            match = nonDNA.exec(seq);
        }
        matches = matches.slice(0, -2); // remove two last characters
        $('#input-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a><span>Invalid DNA sequence at bases ' +
            matches + '.<br>Your ORF should only contain valid DNA bases (AaTtGgCc).</span></div>');
        return false;
    }

    if (seq.length % 3 !== 0) {
        $('#input-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
            'DNA must come in triplets, its length must be a multiple of 3. Current DNA length is: ' +
            seq.length + '</div>');
        return false;
    }
    var triplet;
    for (i = 0; i < seq.length - 3; i = i + 3) {
        triplet = seq.slice(i, i + 3);
        if (triplet === 'TGA' || triplet === 'TAG' || triplet === 'TAA') {
            $('#input-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                '<span>Invalid DNA sequence.<br>Your ORF should not contain internal stop codons (TGA, TAA, TAG).</span></div>');
            return false;
        }
    }
    if (warnings !== '') {
        $('#input-alerts').html('<div class="alert alert-warning"><a class="close" data-dismiss="alert">×</a><span>Warning!<br>' +
            warnings + 'Errors may be caused.</span></div>');
    }
    //The seq string validates maybe with some warnings
    return seq;
}


function draw_coils(ind) {
    $('#sequences-labels').append('<span>' + CCD.gapped_predictions[ind][0] + '</span><br>');
    var coils_pred = CCD.gapped_predictions[ind][1];
    var coils_span = '';
    var i;

    for (i = 0; i < coils_pred.length; i += 1) {
        if (coils_pred[i] === '5') {
            coils_span += '<span style="cursor:default; background-color:rgb(230, 230, 255);">' + coils_pred[i] + '</span>';
        } else if (coils_pred[i] === '6') {
            coils_span += '<span style="cursor:default; background-color:rgb(204, 204, 255);">' + coils_pred[i] + '</span>';
        } else if (coils_pred[i] === '7') {
            coils_span += '<span style="cursor:default; background-color:rgb(179, 179, 255);">' + coils_pred[i] + '</span>';
        } else if (coils_pred[i] === '8') {
            coils_span += '<span style="cursor:default; background-color:rgb(153, 153, 255);">' + coils_pred[i] + '</span>';
        } else if (coils_pred[i] === '9') {
            coils_span += '<span style="cursor:default; background-color:rgb(128, 128, 255);">' + coils_pred[i] + '</span>';
        } else if (coils_pred[i] === '0') {
            coils_span += '<span style="cursor:default; background-color:rgb(43, 43, 255); color: rgb(255, 255, 255)">' + coils_pred[i] + '</span>';
        } else {
            coils_span += '<span style="cursor:default;">' + coils_pred[i] + '</span>';
        }
    }
    $('#sequences-viewer').append('<span style="white-space:pre;">' + coils_span + '</span><br>');
}


function draw_domains(ind) {
    $('#sequences-labels').append('<span>' + CCD.gapped_predictions[ind][0] + '</span><br>');
    var domain_pred = CCD.gapped_predictions[ind][1];
    var domain_span = '';
    var i;

    for (i = 0; i < domain_pred.length; i += 1) {
        if (domain_pred[i] === 'd') {
            domain_span += '<span style="cursor:default; background-color:rgb(244, 212, 222); color: rgb(238, 55, 77);">d</span>';
        } else if (domain_pred[i] === 'G') {
            domain_span += '<span style="cursor:default; background-color:rgb(79, 183, 72); color: rgb(255, 255, 255);">G</span>';
        } else {
            domain_span += '<span style="cursor:default;">' + domain_pred[i] + '</span>';
        }
    }
    $('#sequences-viewer').append('<span style="white-space:pre;">' + domain_span + '</span><br>');
}


function draw_he(ind) {
    $('#sequences-labels').append('<span>' + CCD.gapped_predictions[ind][0] + '</span><br>');
    var he_pred = CCD.gapped_predictions[ind][1];
    var he_span = '';
    var i;

    for (i = 0; i < he_pred.length; i += 1) {
        if (he_pred[i] === 'h') {
            he_span += '<span style="cursor:default; background-color:rgb(44, 171, 212);">h</span>';
        } else if (he_pred[i] === 'e') {
            he_span += '<span style="cursor:default; background-color:rgb(203, 149, 113);">e</span>';
        } else {
            he_span += '<span style="cursor:default;">' + he_pred[i] + '</span>';
        }
    }
    $('#sequences-viewer').append('<span style="white-space:pre;">' + he_span + '</span><br>');
}


function draw_smart(ind) {
    $('#sequences-labels').append('<span>' + CCD.gapped_predictions[ind][0] + '</span><br>');
    var domain_pred = CCD.gapped_predictions[ind][1];
    var transmembrane_str = CCD.gapped_predictions[ind][2];
    var domain_span = '';
    var i;

    for (i = 0; i < domain_pred.length; i += 1) {
        if (transmembrane_str.charAt(i) === 't') {
            domain_span += '<span style="cursor:default; background-color:rgb(255, 165, 0);">' + domain_pred[i] + '</span>';
        } else {
            domain_span += '<span style="cursor:default;">' + domain_pred[i] + '</span>';
        }
    }
    $('#sequences-viewer').append('<span style="white-space:pre;">' + domain_span + '</span><br>');
}


function draw_wait(service_name) {
    $('#sequences-labels').append(service_name + '<br>');
    $('#sequences-viewer').append('Processing...<br>');
}


function raw_predictions() {
    $('#sequences-labels').html('');
    $('#sequences-viewer').html('');
    //draw errors first
    var formatted_errors = '';
    var i;

    if (CCD.pred_errors.length !== 0) {
        for (i = 0; i < CCD.pred_errors.length; i += 1) {
            formatted_errors += CCD.pred_errors[i] + ' error<br>';
        }
        $('#predictions-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a><span>' +
            formatted_errors + '</span></div>');
    }
    $('.service-label').each(function() {
        var prediction_drawn = false;

        if ($(this).children('input').is(':checked') === true) {
            for (i = 0; i < CCD.pred_errors.length; i += 1) {
                if (CCD.pred_errors[i] === $(this).text().substr(3)) { // if there is an error in the name of this prediction name..
                    return true; // continues the .each function to the next element!
                }
            }
            if (CCD.predictions) {
                draw_aa();
                CCD.gapped_predictions = CCD.predictions;
            }
            // loop through the predictions that have arrived and dottified, if there draw it
            for (i = 0; i < CCD.gapped_predictions.length; i += 1) {
                if (CCD.gapped_predictions[i][0] === $(this).text().substr(3)) { // if prediction has arrived
                    switch ($(this).text().substr(3)) {
                        case 'HNN':
                        case 'MLR':
                        case 'DPM':
                        case 'PREDATOR':
                            draw_he(i);
                            break;
                        case 'IUPRED':
                        case 'GLOBPLOT':
                            draw_domains(i);
                            break;
                        case 'SMART':
                            draw_smart(i);
                            break;
                        case 'COILS':
                            draw_coils(i);
                            break;
                        case 'PDB_95':
                        case 'PDB_50to95':
                        case 'PDB_30to50':
                            draw_pdb(i)
                            break;
                    }
                    prediction_drawn = true;
                }
            }
            // if nothing drawn, draw wait since it is checked
            if (prediction_drawn === false) {
                draw_wait($(this).text().substr(3));
            }
        }
    });
}


function arrow_toggle(anchor) {
    var current_arrow = $(anchor).html();
    $(anchor).html(current_arrow === 'more' ? 'less' : 'more'); // php style if
}


function initialize_reset() {
    CCD.starts = [];
    CCD.stops = [];
}

function draw_alignment_table(aligned_isoforms) {
    var formatted_labels = '';
    var formatted_aligned_seqs = '';
    var curr_length = 0;
    var indexes_string = '';
    var alignment_table = '';

    aligned_isoforms.sort(function(a, b) { // sort isoforms based on the name
        var keyA = a[0],
            keyB = b[0];
        // Compare the 2 key strings
        if (keyA < keyB) return -1;
        if (keyA > keyB) return 1;
        return 0;
    });

    while (curr_length < aligned_isoforms[0][1].length - 1) {
        if (curr_length % 10 === 9) {
            indexes_string += '|' + (curr_length + 1).toString();
        } else {
            indexes_string += ' ';
        }
        curr_length = indexes_string.length;
    }
    // Coloring by identity
    var row = '';
    var row_span = '';
    var column = '';
    var all_characters_equality = false;
    var i, j, k;

    aligned_isoforms = aligned_isoforms
    for (i = 0; i < aligned_isoforms.length; i += 1) {
        formatted_labels += aligned_isoforms[i][0] + '<br>';
        row = aligned_isoforms[i][1].replace(/-/g, ' ');
        row_span = '';
        for (j = 0; j < row.length; j += 1) {
            // core of painting
            column = ''; // insert the first character at hand, can this be ommited ?
            for (k = 0; k < aligned_isoforms.length; k += 1) {
                if (aligned_isoforms[k][1][j] !== '-') { // add all characters but gaps
                    column += aligned_isoforms[k][1][j];
                }
            }

            all_characters_equality = column.split('').every(char => char === row.charAt(j)); // every char in the array is equal to the character at hand
            if (all_characters_equality && column.length > 1) {
                row_span += '<span style="cursor:default; background-color: rgb(153, 235, 255);">' + row.charAt(j) + '</span>';
            } else {
                row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
            }
        }
        formatted_aligned_seqs += row_span + '<br>';
    }
    // Coloring end

    alignment_table = '<b>Isoforms alignment</b>&nbsp;&nbsp;<input type="checkbox" id="isoforms-alignment-checkbox" onclick="toggle_isoforms_alignment()"> Show' +
        '<div class="panel-collapse collapse" id="isoforms-alignment-panel"><table class="pred-panel">' +
        '<tr>' +
        '<td class="pred-labels"><br>' + formatted_labels + '</td>' +
        '<td><div class="pred-viewer"><span style="white-space:pre;">' + indexes_string + '</span><br>' +
        '<span style="white-space:pre;">' + formatted_aligned_seqs + '</span></div></td>' +
        '</tr>' +
        '</table></div>';
    return alignment_table;
}
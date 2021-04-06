'use strict';
var CCD = CCD || { // use global CCD as a namespace, window.CCD = window.CCD || {}; would also work. In both cases it is attached to the window object.
    GDPR_accepted: null,
    siteurl: null, // Declaration without assignment is not allowed here, consequently we use null.
    bac_vectors: null,
    ins_vectors: null,
    mam_vectors: null,
    restr_enz: null,
    isoforms_id_dna: null,
    AAseq: null,
    starts: null,
    stops: null,
    vstarts: null,
    vstops: null,
    predictions: null,
    query_str: null,
    gapped_predictions: null,
    pred_errors: null,
    unisearch_ajax: null,
};


function getCookie(cname) {
    var name = cname + '=';
    var ca = document.cookie.split(';');
    for (var i = 0; i < ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0) == ' ') {
            c = c.substring(1);
        }
        if (c.indexOf(name) == 0) {
            return c.substring(name.length, c.length);
        }
    }
    return '';
}


$(document).ready(function() {
    if ($('#first-div').data('title') === 'main_ccd_application') {
        if (getCookie('GDPR_accepted') !== 'yes') {
            CCD.GDPR_accepted = confirm(`General Data Protection Regulation Notice:
As part of the normal operation of CCD, we log the IP address of each CCD user. We use this information only to determine the number and nation of origin of our users, in order to get additional funding for development of this and other services. CCD does not collect or store any other type of information about you. By clicking "OK", you agree to the above. Your preference will be stored in a cookie on your browser.`);
            if (CCD.GDPR_accepted) {
                document.cookie = "GDPR_accepted=yes; expires=Thu, 18 Dec 3000 12:00:00 UTC";
            } else {
                document.write('We are sorry but you cannot use our service.');
                return '';
            }
        }
        CCD.siteurl = window.location.href;
        $('[data-toggle="tooltip"]').tooltip({
            html: true,
            trigger: 'hover'
        });
        select_all_serv();
        // Ajax requests abort facilitate
        jQuery.xhrPool = [];
        jQuery.xhrPool.abortAll = function() {
            var requests = [];
            var index;

            for (index in this) {
                if (isFinite(index) === true) {
                    requests.push(this[index]);
                }
            }
            for (index in requests) {
                if (isFinite(index) === true) {
                    requests[index].abort();
                }
            }
        };
        $(document).ajaxSend(function(event, jqXHR, options) {
            jQuery.xhrPool.push(jqXHR);
        });
        $(document).ajaxComplete(function(event, jqXHR, options) {
            jQuery.xhrPool.remove(jqXHR);
        });
        jQuery.xhrPool.remove = function(jqXHR) {
            var index;

            for (index in this) {
                if (this[index] === jqXHR) {
                    jQuery.xhrPool.splice(index, 1);
                    break;
                }
            }
        };
        hide_byclass('section1');

        $('.collapse.in').prev('.panel-heading').addClass('active');
        $('#bs-collapse')
            .on('show.bs.collapse', function(a) {
                $(a.target).prev('.panel-heading').addClass('active');
            })
            .on('hide.bs.collapse', function(a) {
                $(a.target).prev('.panel-heading').removeClass('active');
            });
        $('#search-id').on('keyup', function(e) {
            if (e.keyCode == 13) {
                search_uniprot();
            }
        });
    }
});


function search_uniprot() {
    var uniprot_id = $('#search-id').val();
    uniprot_id = uniprot_id.trim(); // whitespace removal
    clear_uniprot_search();

    if (validate_search_uniprot_input(uniprot_id) === false) {
        return false;
    }

    clear_all();
    hide_byclass('section1');
    jQuery.xhrPool.abortAll();

    // initialize and reset
    initialize_reset();

    $('.uniprot-search-progress-bar').css('display', 'inline');
    CCD.unisearch_ajax = $.ajax({
        type: 'GET',
        url: CCD.siteurl + 'search_uniprot/' + uniprot_id,
        async: true,
        success: function(response) {
            clear_uniprot_search(); // used to hide progress bar
            if (response.status === 'ERROR') {
                $('#uniprot-id-search-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                    '<span>' + response.reason + '</span></div>');
                return false;
            }

            var protein_isoforms_list = response.protein_isoforms_list; // no need to parse cause it is already json
            var formatted_search_results = '';
            var isoform_description = '';
            var isoform_matches = '';
            var cds_link = '';
            var entry_link = '';
            var delimiter;
            var text_color = '';
            var mismatch_note = '';
            var disable = '';
            var alignment_table = '';
            var i, j, k;

            CCD.isoforms_id_dna = [];

            protein_isoforms_list.sort(function(a, b) {
                var keyA = a.isoform_id,
                    keyB = b.isoform_id;
                // Compare the 2 key strings
                if (keyA < keyB) return -1;
                if (keyA > keyB) return 1;
                return 0;
            });

            if (response.aligned_isoforms) {
                alignment_table = draw_alignment_table(response.aligned_isoforms);
            }

            for (i = 0; i < protein_isoforms_list.length; i += 1) {
                isoform_description = '';
                delimiter = '';
                text_color = '';
                mismatch_note = '';
                disable = '';

                for (j = 0; j < protein_isoforms_list[i].isoform_description.length - 1; j += 1) {
                    isoform_description += delimiter + protein_isoforms_list[i].isoform_description[j + 1].slice(0, -1);
                    delimiter = '<br>';
                }
                if (protein_isoforms_list[i].isoform_mismatches !== '0') {
                    text_color = 'color: orange;';
                    mismatch_note = 'This isoform does not have an exact match.';
                }
                if (protein_isoforms_list[i].isoform_dna === 'None') {
                    text_color = 'color: red;';
                    disable = 'disabled';
                }
                formatted_search_results += '<span style="display: inline-block; width: 630px;"><label class="label-poly-check">' +
                    '<input type="radio" name="isoform-radio" class="isoform-check" onclick="isoform_dna()" ' + disable + '> <b style="' + text_color + '">' +
                    protein_isoforms_list[i].isoform_id + '</b></label>' +
                    'Length: ' + protein_isoforms_list[i].isoform_protein_length + '&nbsp;&nbsp;' +
                    'Short protein sequence: <span style=" font-family: monospace;">' + protein_isoforms_list[i].isoform_short_protein + '</span></span>&nbsp;&nbsp;' +
                    '<a data-toggle="collapse" data-target="#isoform' + (i + 1) + '-info" onclick="arrow_toggle(this);">more</a><br>' +
                    '<div id="isoform' + (i + 1) + '-info" class="panel-collapse collapse">' + protein_isoforms_list[i].isoform_description[0] // first line inindented!
                    +
                    '<div class="extra-margin-left">' + isoform_description + '</div>';

                if (protein_isoforms_list[i].isoform_matches !== '0') {
                    isoform_matches = '';
                    delimiter = '';

                    for (k = 0; k < protein_isoforms_list[i].isoform_matched_with.length; k += 1) {
                        cds_link = '<a href="' + protein_isoforms_list[i].isoform_matched_with[k].cds_url +
                            '" target="_blank">' + protein_isoforms_list[i].isoform_matched_with[k].cds_id + '</a>';
                        isoform_matches += delimiter + cds_link;
                        delimiter = ', '; // Each time except the first we will have the delimiter. Programming pattern:)
                    }
                    formatted_search_results += 'Matches: ' + isoform_matches + '<br>';
                }
                if (protein_isoforms_list[i].isoform_mismatches !== '0') {
                    var subs = '';
                    delimiter = '';
                    for (k = 0; k < protein_isoforms_list[i].isoform_mismatched_with.length; k += 1) {
                        subs += delimiter + protein_isoforms_list[i].isoform_mismatched_with[k][0] + '->' + protein_isoforms_list[i].isoform_mismatched_with[k][1];
                        delimiter = ',&nbsp;&nbsp;';
                    }
                    formatted_search_results += 'Substitutions: ' + subs + '<br>';
                }
                formatted_search_results += '<br></div>';

                CCD.isoforms_id_dna.push([protein_isoforms_list[i].isoform_id, protein_isoforms_list[i].isoform_dna]);
            }


            entry_link = '<a href="http://www.uniprot.org/uniprot/' + response.protein_id +
                '" id="uniprot-id">' + response.protein_id + '</a>';
            formatted_search_results = '<br><div class="alert alert-info"><a class="close" data-dismiss="alert">×</a><span><strong>Tip:</strong> ' +
                'Once you click on an isoform, its DNA will be automatically pasted on the DNA sequence area below.</span></div>' +
                alignment_table +
                '<div style="border-bottom-style: solid; border-bottom-width: 1px;"><br><b>Uniprot ID | ' +
                entry_link + '</b><br></div><br>Organism: <i>' + response.protein_source_organism +
                ' (' + response.protein_source_organism_common + ')</i>' +
                '<br>This entry describes the following isoforms:<br>' +
                formatted_search_results + '</div>';
            $('#uniprot-id-search-results').html(formatted_search_results);
        },
        error: function(jqXHR, exception) {
            clear_uniprot_search(); // used to hide progress bar
            if (exception === 'abort') {

            } else {
                $('#uniprot-id-search-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                    '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
            }
        }
    });
}


function isoform_dna() {
    var i;

    for (i = 0; i < CCD.isoforms_id_dna.length; i += 1) {
        if (CCD.isoforms_id_dna[i][0] === $('input[name="isoform-radio"]:checked').next().text()) {
            document.getElementById('dna-orf').value = CCD.isoforms_id_dna[i][1];
        }
    }
}


function get_predictions() {
    $('#aa-viewer').css('overflow', 'auto');
    var DNAseq = document.getElementById('dna-orf').value;

    DNAseq = validateDNA(DNAseq);
    if (DNAseq === false) {
        return false;
    }

    clear_all();
    hide_byclass('section1');
    show_byclass('section1');
    jQuery.xhrPool.abortAll();

    // initialize and reset
    initialize_reset();

    $.ajax({
        async: false,
        type: 'POST',
        url: CCD.siteurl + 'translate',
        data: {
            DNAseq: DNAseq
        },
        success: function(response) {
            var json_obj = $.parseJSON(response);

            if (json_obj.status === 'VALIDATE_ERROR') {
                $('#input-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                    '<span>Server side validation error.</span></div>');
            } else {
                CCD.AAseq = json_obj.AAseq;
                draw_aa();
            }
        },
        error: function() {
            $('#predictions-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
        }
    });


    $('#predictions-info-placeholder').html('<div class="alert alert-info"><a class="close" data-dismiss="alert">×</a><span><strong>Tips:</strong> ' +
        '1) Click again on a residue to unmark it. ' +
        '2) The same residue can be marked as start and stop. ' +
        '3) Place the cursor over a residue to see its number in the amino acid sequence. ');
    prediction_services();
}


function minus_one(this_obj) {
    $('.' + $(this_obj).parent().attr('class')).css('display', 'none');
}


function reset_malignment() {
    $('[class^="malignment_"]').css('display', 'inline');
}


function prediction_services() {
    //Runs the predictions that are listed in the prediciton servers button (index.html), if checked
    //The prediction are named exactly as stated in the label
    //The backend returns array of [prediction_name, results_string] tuples
    $('#predictions-alerts-placeholder').html('');
    CCD.predictions = [];
    CCD.pred_errors = [];

    draw_predictions();
    $('.service-label').each(function() {
        if ($(this).children('input').is(':checked') === true) {
            $.ajax({
                async: true,
                type: 'GET',
                url: CCD.siteurl + 'prediction/' + $(this).text().substr(3),
                success: function(response) {
                    var json_obj = $.parseJSON(response);
                    var prediction = json_obj.prediction;
                    var transmembrane_str = '';
                    var flag2continue = false; // to avoid using loopnext: label, which is the JS version of evil GOTO
                    var i, j;
                    if (json_obj.error === 'NO') {
                        if (prediction[0] === 'SMART') {
                            for (i = 0; i < prediction[1].length; i += 1) {
                                flag2continue = false;
                                for (j = 0; j < prediction[2].length; j += 1) {
                                    if (i >= prediction[2][j][0] - 1 && i < prediction[2][j][1] && prediction[2][j][2] === 'transmembrane region') { // we only care to display transmembrane regions
                                        transmembrane_str += 't';
                                        flag2continue = true;
                                    }
                                }
                                if (flag2continue) {
                                    continue;
                                } else {
                                    transmembrane_str += '-';
                                }
                            }
                            prediction[2] = transmembrane_str;
                        }
                        CCD.predictions.push(prediction);
                        $('#aa-viewer').css('overflow', 'hidden'); // to hide the scrool bar
                    }
                    if (json_obj.error === 'YES') {
                        CCD.pred_errors.push(prediction);
                    }
                    draw_predictions();
                },
                error: function(jqXHR, exception) {
                    if (exception === 'abort') {
                        $('#predictions-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                            '<span>Predictions request aborted.</span></div>');
                    } else {
                        $('#predictions-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                            '<span>Communication problem with CCD 2 server.<br>' +
                            'Please try again or reload the page.</span></div>');
                    }
                }
            });
        }
    });
}


function onclickAA(ev) {
    var current_state = $('#mark-radios label.active input').val();

    // coloring and stroring data
    switch (ev.target.style.backgroundColor) {
        case 'lawngreen': // already marked as a start
            if (current_state === 'starts-mode') { // if markmode starts
                ev.target.style.backgroundColor = ''; // deselect
                CCD.starts.splice(CCD.starts.indexOf(ev.target.getAttribute('name')), 1); // remove from starts array
            }
            if (current_state === 'stops-mode') { // if markmode stops
                ev.target.style.backgroundColor = 'yellow'; // make it both start and stop
                CCD.stops.push(ev.target.getAttribute('name')); // add to stops array
            }
            break;
        case 'tomato': // already marked as a stop
            if (current_state === 'starts-mode') { // if markmode starts
                ev.target.style.backgroundColor = 'yellow'; // make it both start and stop
                CCD.starts.push(ev.target.getAttribute('name')); // add to starts array
            }
            if (current_state === 'stops-mode') { // if markmode stops
                ev.target.style.backgroundColor = ''; // deselect
                CCD.stops.splice(CCD.stops.indexOf(ev.target.getAttribute('name')), 1); // remove from stops array
            }
            break;
        case 'yellow': // already marked as a start and a stop
            if (current_state === 'starts-mode') { // if markmode starts
                ev.target.style.backgroundColor = 'tomato'; // select only as stop
                CCD.starts.splice(CCD.starts.indexOf(ev.target.getAttribute('name')), 1); // remove from starts array
            }
            if (current_state === 'stops-mode') { // if markmode stops
                ev.target.style.backgroundColor = 'lawngreen'; // select only as start
                CCD.stops.splice(CCD.stops.indexOf(ev.target.getAttribute('name')), 1); // remove from stops array
            }
            break;
        case '': // if no marking yet
            if (current_state === 'starts-mode') { // if markmode starts
                ev.target.style.backgroundColor = 'lawngreen'; // select as a start
                CCD.starts.push(ev.target.getAttribute('name')); // add to starts array
            }
            if (current_state === 'stops-mode') { // if markmode stops
                ev.target.style.backgroundColor = 'tomato'; // select as a stop
                CCD.stops.push(ev.target.getAttribute('name')); // add to stops array
            }
            break;
    }
}


function save_predictions() {
    var link = document.createElement('a');
    var mimeType = 'text/plain';
    var pfname = $('#predictions-filename').val();
    var formatted_predictions = 'AAsequence,' + CCD.AAseq + '\n';
    var i;

    if (pfname === null || pfname === '') {
        pfname = 'ccd_predictions.csv';
    } else {
        pfname += '.csv';
    }
    link.setAttribute('download', pfname);
    for (i = 0; i < CCD.predictions.length; i += 1) {
        formatted_predictions += CCD.predictions[i][0] + ',';
        formatted_predictions += CCD.predictions[i][1] + '\n';
    }
    link.setAttribute('href', 'data:' + mimeType + ';charset=utf-8,' + encodeURIComponent(formatted_predictions));
    document.body.appendChild(link); // Firefox requires the link to be in the body
    link.click();
    document.body.removeChild(link); // remove the link when done
}


function seq_file() {
    var f = document.getElementById('file-seq').files[0];

    if (f) {
        var r = new FileReader();
        r.onload = function(e) {
            var contents = e.target.result;
            document.getElementById('dna-orf').value = contents;
        };
        r.readAsText(f);
    }
}


function sample_sequence() {
    var dna = '>ENA|AAH05389|AAH05389.1 Homo sapiens (human) geminin, DNA replication inhibitor\n' +
        'ATGAATCCCAGTATGAAGCAGAAACAAGAAGAAATCAAAGAGAATATAAAGAATAGTTCT' +
        'GTCCCAAGAAGAACTCTGAAGATGATTCAGCCTTCTGCATCTGGATCTCTTGTTGGAAGA' +
        'GAAAATGAGCTGTCCGCAGGCTTGTCCAAAAGGAAACATCGGAATGACCACTTAACATCT' +
        'ACAACTTCCAGCCCTGGGGTTATTGTCCCAGAATCTAGTGAAAATAAAAATCTTGGAGGA' +
        'GTCACCCAGGAGTCATTTGATCTTATGATTAAAGAAAATCCATCCTCTCAGTATTGGAAG' +
        'GAAGTGGCAGAAAAACGGAGAAAGGCGCTGTATGAAGCACTTAAGGAAAATGAGAAACTT' +
        'CATAAAGAAATTGAACAAAAGGACAATGAAATTGCCCGCCTGAAAAAGGAGAATAAAGAA' +
        'CTGGCAGAAGTAGCAGAACATGTACAGTATATGGCAGAGCTAATAGAGAGACTGAATGGT' +
        'GAACCTCTGGATAATTTTGAATCACTGGATAATCAGGAATTTGATTCTGAAGAAGAAACT' +
        'GTTGAGGATTCTCTAGTGGAAGACTCAGAAATTGGCACGTGTGCTGAAGGAACTGTATCT' +
        'TCCTCTACGGATGCAAAGCCATGTATATGA';
    document.getElementById('dna-orf').value = dna;
}


function reset_aa() {
    $('.amino-acid').css('background-color', '');
    CCD.starts = [];
    CCD.stops = [];

}


function select_all_serv() {
    $('.service-label > input').each(function() {
        this.checked = true;
    });
}


function deselect_all_serv() {
    $('.service-label > input').each(function() {
        this.checked = false;
    });
}

function clear_sequence() {
    document.getElementById('dna-orf').value = '';
    document.getElementById('file-seq').type = '';
    document.getElementById('file-seq').type = 'file';
    $('#input-alerts').html('');
    $('.isoform-check').each(function(i) {
        this.checked = false;
    });
}


function save_alignment() {
    //preparing the alignment in fasta format
    var fasta = []
    var filename = ''
    var uniprot_id = $('#search-id').val();
    if (uniprot_id) {
        filename = `${uniprot_id}_alignment.fasta`
    } else {
        filename = 'CCD_alignments.fasta'
    }
    for (var i = 0; i < CCD.multiple_alignment.length; i += 1) {
        fasta.push(`>${CCD.multiple_alignment[i][0]}\n${CCD.multiple_alignment[i][1]}\n`);
    }
    var file = new Blob(fasta, { type: 'text/plain' });
    var a = document.createElement("a"),
        url = URL.createObjectURL(file);
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    setTimeout(function() {
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);
    }, 0);
}


function toggle_isoforms_alignment() {
    var alignment_visibility = $('#isoforms-alignment-checkbox').is(':checked');
    if (alignment_visibility === true) {
        // make area visible
        $('#isoforms-alignment-panel').css('display', 'table');
    } else if (alignment_visibility === false) {
        $('#isoforms-alignment-panel').css('display', 'none');
    } else {
        alert('Toggling multiple alignment failed.');
    }
}
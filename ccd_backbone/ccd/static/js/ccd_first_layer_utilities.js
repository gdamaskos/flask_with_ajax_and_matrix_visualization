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
    polypeptides: null,
    primers: null,
    tagged_cleaved: null,
    proteases: null,
    tag_ajax: null,
    cleave_ajax: null,
    crysol_pending: null,
    unisearch_ajax: null,
    multiple_alignment: null,
    crysol_raw: null,
    crysol_tagged: null,
    crysol_cleaved: null,
    crysol_raw_sort_switch: null,
    crysol_tagged_sort_switch: null,
    crysol_cleaved_sort_switch: null,
    crysol_task_ids: null,
    crysol_stop_recursion: null
};


window.onbeforeunload = function (e) { // revoke any celery tasks that are pending 
    if (CCD.crysol_task_ids) {
        crysol_revoke();
    }
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


$(document).ready(function () {
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
        $.ajax({
            async: false,
            url: CCD.siteurl + 'get_vectors',
            type: 'GET',
            success: function (response) {
                var i;
                // these variables are global
                CCD.bac_vectors = response.bacterial_vectors;
                CCD.ins_vectors = response.insect_vectors;
                CCD.mam_vectors = response.mammalian_vectors;
                CCD.restr_enz = response.restriction_enzymes;

                $('#vectors-list').append('<option disabled="true">Bacterial</option>');
                for (i = 0; i < response.bacterial_vectors.length; i += 1) {
                    $('#vectors-list').append('<option>' + response.bacterial_vectors[i].name + '</option>');
                }
                $('#vectors-list').append('<option disabled="true">Insect</option>');
                for (i = 0; i < response.insect_vectors.length; i += 1) {
                    $('#vectors-list').append('<option>' + response.insect_vectors[i].name + '</option>');
                }
                $('#vectors-list').append('<option disabled="true">Mammalian</option>');
                for (i = 0; i < response.mammalian_vectors.length; i += 1) {
                    $('#vectors-list').append('<option>' + response.mammalian_vectors[i].name + '</option>');
                }
                for (i = 0; i < response.restriction_enzymes.length; i += 1) {
                    $('#restr-fw-overhang').append('<option>' + response.restriction_enzymes[i].name + '</option>');
                    $('#restr-rv-overhang').append('<option>' + response.restriction_enzymes[i].name + '</option>');
                }
                $('#restr-fw-overhang').append('<option>' + 'Other..' + '</option>');
                $('#restr-rv-overhang').append('<option>' + 'Other..' + '</option>');
            },
            error: function () {
                alert('Communication problem with CCD 2 server.\nPlease reload the page.');
            }
        });

        $('#vectors-list').bind('keyup', function () {
            $(this).trigger('change');
        });
        enable_lic();
        select_all_serv();
        // scroll synchronization of aa-viewer and evolutionator-viewer, their bars will be hidden on success
        $('#sequences-viewer').on('scroll', function () {
            $('#aa-viewer').scrollLeft($(this).scrollLeft());
            $('#evolutionator-viewer').scrollLeft($(this).scrollLeft());
        });
        $('#evolutionator-viewer').on('scroll', function () {
            $('#sequences-viewer').scrollLeft($(this).scrollLeft());
        });
        // horizontal scroll synchronization for label and sequences in multiple alignment panels
        $('#evolutionator-labels').on('scroll', function () {
            $('#evolutionator-viewer').scrollTop($(this).scrollTop());
        });
        $('#evolutionator-viewer').on('scroll', function () {
            $('#evolutionator-labels').scrollTop($(this).scrollTop());
        });
        $('#evolutionator-labels').css('overflow-y', 'hidden'); // to hide the scrool bar
        // Ajax requests abort facilitate
        jQuery.xhrPool = [];
        jQuery.xhrPool.abortAll = function () {
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
        $(document).ajaxSend(function (event, jqXHR, options) {
            jQuery.xhrPool.push(jqXHR);
        });
        $(document).ajaxComplete(function (event, jqXHR, options) {
            jQuery.xhrPool.remove(jqXHR);
        });
        jQuery.xhrPool.remove = function (jqXHR) {
            var index;

            for (index in this) {
                if (this[index] === jqXHR) {
                    jQuery.xhrPool.splice(index, 1);
                    break;
                }
            }
        };
        hide_byclass('section1');
        hide_byclass('section2');
        hide_byclass('section3');
        hide_byclass('tag');
        hide_byclass('cleave');
        $('.collapse.in').prev('.panel-heading').addClass('active');
        $('#bs-collapse')
            .on('show.bs.collapse', function (a) {
                $(a.target).prev('.panel-heading').addClass('active');
            })
            .on('hide.bs.collapse', function (a) {
                $(a.target).prev('.panel-heading').removeClass('active');
            });
        $('#search-id').on('keyup', function (e) {
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
    hide_byclass('section2');
    hide_byclass('section3');
    hide_byclass('tag');
    hide_byclass('cleave');
    jQuery.xhrPool.abortAll();

    // initialize and reset
    initialize_reset();

    $('.uniprot-search-progress-bar').css('display', 'inline');
    CCD.unisearch_ajax = $.ajax({
        type: 'GET',
        url: CCD.siteurl + 'search_uniprot/' + uniprot_id,
        async: true,
        success: function (response) {
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
            var isoform_mismatches = '';
            var cds_link = '';
            var entry_link = '';
            var delimiter;
            var text_color = '';
            var mismatch_note = '';
            var disable = '';
            var alignment_table = '';
            var i, j, k;

            CCD.isoforms_id_dna = [];

            protein_isoforms_list.sort(function (a, b) {
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
        error: function (jqXHR, exception) {
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
    hide_byclass('section2');
    hide_byclass('section3');
    hide_byclass('tag');
    hide_byclass('cleave');
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
        success: function (response) {
            var json_obj = $.parseJSON(response);

            if (json_obj.status === 'VALIDATE_ERROR') {
                $('#input-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                    '<span>Server side validation error.</span></div>');
            } else {
                CCD.AAseq = json_obj.AAseq;
                draw_aa();
            }
        },
        error: function () {
            $('#predictions-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
        }
    });

    evolutionator_call();

    $('#predictions-info-placeholder').html('<div class="alert alert-info"><a class="close" data-dismiss="alert">×</a><span><strong>Tips:</strong> ' +
        '1) Click again on a residue to unmark it. ' +
        '2) The same residue can be marked as start and stop. ' +
        '3) Place the cursor over a residue to see its number in the amino acid sequence. ' +
        '4) Click on the &#10134; symbol before a multiple alignment label to remove it from the view, press Reset in Alignment panel to redraw all.</span></div>');
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
    $('.service-label').each(function () {
        if ($(this).children('input').is(':checked') === true) {
            $.ajax({
                async: true,
                type: 'GET',
                url: CCD.siteurl + 'prediction/' + $(this).text().substr(3),
                success: function (response) {
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
                error: function (jqXHR, exception) {
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


function get_primers_peptides() {
    var starts = CCD.starts;
    var stops = CCD.stops;
    var method_value = [];

    clear_after_primers();
    hide_byclass('tag');
    unpress_tag_button();
    hide_byclass('cleave');
    unpress_cleave_button();
    hide_byclass('section3');
    $('#poly-check-select-pp').prop('checked', false); // Uncheck the Select all/none box 
    // Uncheck all the boxes
    $('.poly-check1').each(function () {
        this.checked = false;
    });
    $('.poly-check2').each(function () {
        this.checked = false;
    });
    $('.poly-check3').each(function () {
        this.checked = false;
    });


    if (starts.length === 0 && stops.length === 0) {
        $('#methods-overhangs-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
            '<span>No start or stop is marked.</span></div>');
        return false;
    }

    show_byclass('section2');

    if (document.getElementById('annealing_temperature_radio').checked === true) {
        method_value = [0, document.getElementById('annealing_temperature_holder').value];
    } else if (document.getElementById('base_length_radio').checked === true) {
        method_value = [1, document.getElementById('base_length_holder').value];
    }

    $.ajax({
        async: false,
        url: CCD.siteurl + 'primersnpeptides',
        data: {
            starts: JSON.stringify(starts),
            stops: JSON.stringify(stops),
            method_value: JSON.stringify(method_value),
            FWoverhang: $('#fw-overhang').text(),
            RVoverhang: $('#rv-overhang').text()
        },
        type: 'POST',
        success: function (response) {
            var json_obj = $.parseJSON(response);
            var primers_viewer = $('#primers-viewer');
            var polypeptides_labels = $('#polypeptides-labels');
            var polypeptides_viewer = $('#polypeptides-viewer');
            var MWs_viewer = $('#MWs-viewer');
            var pIs_viewer = $('#pIs-viewer');
            var epsilons_viewer = $('#epsilons-viewer');
            CCD.vstarts = json_obj.vstarts;
            CCD.vstops = json_obj.vstops;
            var MWs = json_obj.MWs;
            var pIs = json_obj.pIs;
            var epsilons = json_obj.epsilons;
            var i;
            // so that a global variable will be available for save
            CCD.primers = json_obj.primers;
            build_polypeptides(); // this creates a global variable with all polypeptides in it
            var primers_errors = '';


            polypeptides_labels.html('');
            polypeptides_viewer.html('');
            MWs_viewer.html('');
            pIs_viewer.html('');
            epsilons_viewer.html('');
            primers_viewer.html('');

            if (json_obj.status === 'VALIDATE_ERROR') {
                $('#methods-overhangs-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                    '<span>Server side validation error. An overhang can only contain letters.</span></div>');
                return false;
            }
            for (i = 0; i < CCD.primers.length; i += 1) {
                primers_viewer.append('<tr><td style="border: 1px solid; font-weight:normal; font-family:serif; font-size:17px; overflow:hidden; white-space:nowrap; width:16%;">' +
                    $('#pp-name').val() + CCD.primers[i][0] +
                    '</td><td style="border: 1px solid; font-family:monospace; font-size:17px; white-space: nowrap;"><div style="overflow:auto; width:100%">' +
                    CCD.primers[i][1] + '</div></td></tr>');
            }
            $('#primers-alerts-placeholder').html('');
            for (i = 0; i < json_obj.errors.length; i += 1) {
                primers_errors += json_obj.errors[i] + '<br>';
            }
            if (json_obj.errors.length !== 0) {
                primers_errors = '';
                for (i = 0; i < json_obj.errors.length; i += 1) {
                    primers_errors += json_obj.errors[i] + '<br>';
                }
                $('#primers-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a><span>' +
                    primers_errors + '</span></div>');
            }
            for (i = 0; i < CCD.polypeptides.length; i += 1) {
                polypeptides_labels.append('<label class="label-poly-check"><input type="checkbox" class="poly-check1" id="">' +
                    CCD.polypeptides[i][0] + '</label><br>');
                MWs_viewer.append(MWs[i] + '<br>');
                pIs_viewer.append(pIs[i] + '<br>');
                epsilons_viewer.append(epsilons[i] + '<br>');
                polypeptides_viewer.append('<span style="white-space:pre;">' + CCD.polypeptides[i][1] + '</span><br>');
            }
        },
        error: function () {
            $('#primers-alerts-placeholder').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
        }
    });
    plasmid_populate();
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


function save_primers() {
    var link = document.createElement('a');
    var mimeType = 'text/plain';
    var pfname = $('#primers-filename').val();
    var formatted_primers = '';
    var i;

    if (pfname === null || pfname === '') {
        pfname = 'ccd_primers.csv';
    } else {
        pfname += '.csv';
    }
    link.setAttribute('download', pfname);
    for (i = 0; i < CCD.primers.length; i += 1) {
        formatted_primers += CCD.primers[i][0] + ',';
        formatted_primers += CCD.primers[i][1] + '\n';
    }
    link.setAttribute('href', 'data:' + mimeType + ';charset=utf-8,' + encodeURIComponent(formatted_primers));
    document.body.appendChild(link); // Firefox requires the link to be in the body
    link.click();
    document.body.removeChild(link); // remove the link when done
}


function save_polypeptides() {
    var link = document.createElement('a');
    var mimeType = 'text/plain';
    var pfname = $('#polypeptides-filename').val();
    var formatted_polypeptides = '';
    var i, j;

    if (pfname === null || pfname === '') {
        pfname = 'ccd_polypeptides.csv';
    } else {
        pfname += '.csv';
    }
    link.setAttribute('download', pfname);
    for (i = 0; i < CCD.polypeptides.length; i += 1) {
        formatted_polypeptides += CCD.polypeptides[i][0] + ',';
        formatted_polypeptides += CCD.polypeptides[i][1] + '\n';
    }
    if ($('#tagged-pp-title').css('display') === 'block') {
        for (i = 0; i < CCD.tagged_cleaved.length; i += 1) {
            formatted_polypeptides += CCD.tagged_cleaved[i][0] + ' (tagged),';
            formatted_polypeptides += CCD.tagged_cleaved[i][4] + '\n';
        }
    }
    if ($('#cleaved-pp-title').css('display') === 'block') {
        for (i = 0; i < CCD.proteases.length; i += 1) {
            for (j = 0; j < CCD.tagged_cleaved.length; j += 1) {
                formatted_polypeptides += CCD.tagged_cleaved[j][0] + ' (cleaved with ' + CCD.proteases[i] + ' ),';
                formatted_polypeptides += CCD.tagged_cleaved[i][8][i] + '\n';
            }
        }
    }
    link.setAttribute('href', 'data:' + mimeType + ';charset=utf-8,' + encodeURIComponent(formatted_polypeptides));
    document.body.appendChild(link); // Firefox requires the link to be in the body
    link.click();
    document.body.removeChild(link); // remove the link when done
}


function seq_file() {
    var f = document.getElementById('file-seq').files[0];

    if (f) {
        var r = new FileReader();
        r.onload = function (e) {
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


function enable_lic() {
    $('#lic').css('display', 'block');
    $('#rc').css('display', 'none');
    $('#custom').css('display', 'none');
    vectors_update();
}

function enable_rc() {
    $('#lic').css('display', 'none');
    $('#rc').css('display', 'block');
    $('#custom').css('display', 'none');
    restr_update_fw();
    restr_update_fwfiller();
    restr_update_start_codon();
    restr_update_rv();
    restr_update_rvfiller();
    restr_update_stop_codon();
}


function enable_custom() {
    $('#lic').css('display', 'none');
    $('#rc').css('display', 'none');
    $('#custom').css('display', 'block');
    custom_update();
}


function vectors_update() {
    var max_len = Math.max(CCD.bac_vectors.length, CCD.ins_vectors.length, CCD.mam_vectors.length);
    var fw_overhang = '';
    var rv_overhang = '';
    var i;

    for (i = 0; i < max_len; i += 1) {
        if (i < CCD.bac_vectors.length && $('#vectors-list option:selected').text() === CCD.bac_vectors[i].name) {
            fw_overhang = CCD.bac_vectors[i].fw_overhang;
            rv_overhang = CCD.bac_vectors[i].rv_overhang;
            break;
        }
        if (i < CCD.ins_vectors.length && $('#vectors-list option:selected').text() === CCD.ins_vectors[i].name) {
            fw_overhang = CCD.ins_vectors[i].fw_overhang;
            rv_overhang = CCD.ins_vectors[i].rv_overhang;
            break;
        }
        if (i < CCD.mam_vectors.length && $('#vectors-list option:selected').text() === CCD.mam_vectors[i].name) {
            fw_overhang = CCD.mam_vectors[i].fw_overhang;
            rv_overhang = CCD.mam_vectors[i].rv_overhang;
            break;
        }
    }

    $('#fw-overhang').html(fw_overhang);
    $('#rv-overhang').html(rv_overhang);
    mini_validator_fw(fw_overhang);
    mini_validator_rv(rv_overhang);
    end_notification_fw(fw_overhang);
    end_notification_rv(rv_overhang);
    match_vectors();
}


function custom_update() {
    $('#fw-overhang').html($('#custom-fw-overhang').val());
    $('#rv-overhang').html($('#custom-rv-overhang').val());
    mini_validator_fw($('#custom-fw-overhang').val());
    mini_validator_rv($('#custom-rv-overhang').val());
    end_notification_fw($('#custom-fw-overhang').val());
    end_notification_rv($('#custom-rv-overhang').val());
}


function toggle_vectors() {
    var vectors_visibility = $('#vectors-toggler').is(':checked');
    if (vectors_visibility === true) {
        // make matching vectors visible
        $('#matching-vectors').css('display', 'block');
    } else if (vectors_visibility === false) {
        $('#matching-vectors').css('display', 'none');
    } else {
        alert('Toggling matching vectors failed.');
    }
}


function toggle_alignment() {
    var alignment_visibility = $('#alignment-toggler').is(':checked');
    if (alignment_visibility === true) {
        // make area visible
        $('#alignment-panel').css('display', 'table');
        draw_predictions();
    } else if (alignment_visibility === false) {
        $('#alignment-panel').css('display', 'none');
        draw_predictions();
    } else {
        alert('Toggling multiple alignment failed.');
    }
}


function toggle_color() {
    var color_visibility = $('#color-toggler').is(':checked');
    if (color_visibility === true) {
        // draw alignment with color
        draw_malignment();
    } else if (color_visibility === false) {
        // draw alignment without color
        raw_malignment();
    } else {
        alert('Toggling multiple alignment failed.');
    }
}


function restr_update_fw() {
    var index = document.getElementById('restr-fw-overhang').selectedIndex;
    $('#fwcustom-restr-site').css('display', 'none');
    if (index === 13) { // why is that ?
        $('#fwcustom-restr-site').css('display', 'block');
        restr_update_fwcustom();
    } else {
        var fw_overhang = $('#filler-fw-overhang').val() +
            CCD.restr_enz[index].sequence +
            $('#scodon-fw-overhang').val();
        $('#fw-overhang').html(fw_overhang);
        mini_validator_fw(fw_overhang);
        end_notification_fw(fw_overhang);
    }
}


function restr_update_fwfiller() {
    var index = document.getElementById('restr-fw-overhang').selectedIndex;
    if (index === 13) {
        restr_update_fwcustom();
    } // if the selection is to add custom restriction enzyme
    else {
        restr_update_fw();
    }
}


function restr_update_start_codon() {
    var index = document.getElementById('restr-fw-overhang').selectedIndex;
    if (index === 13) {
        restr_update_fwcustom();
    } // if the selection is to add custom restriction enzyme
    else {
        restr_update_fw();
    }
}


function restr_update_fwcustom() {
    var fw_overhang = $('#filler-fw-overhang').val() +
        $('#fwcustom-restr-site').val() +
        $('#scodon-fw-overhang').val();
    var rv_overhang = $('#rv-overhang').html();
    $('#fw-overhang').html(fw_overhang);
    mini_validator_fw(fw_overhang);
    end_notification_fw(fw_overhang);
}


function restr_update_rv() {
    var index = document.getElementById('restr-rv-overhang').selectedIndex;
    $('#rvcustom-restr-site').css('display', 'none');
    if (index === 13) { // if the selection is to add custom restriction enzyme
        $('#rvcustom-restr-site').css('display', 'block');
        restr_update_rvcustom();
    } else {
        var rv_overhang = $('#filler-rv-overhang').val() +
            CCD.restr_enz[index].sequence +
            $('#scodon-rv-overhang').val();
        $('#rv-overhang').html(rv_overhang);
        mini_validator_rv(rv_overhang);
        end_notification_rv(rv_overhang);
    }
}


function restr_update_rvfiller() {
    var index = document.getElementById('restr-rv-overhang').selectedIndex;
    if (index === 13) {
        restr_update_rvcustom();
    } // if the selection is to add custom restriction enzyme
    else {
        restr_update_rv();
    }
}


function restr_update_stop_codon() {
    var index = document.getElementById('restr-rv-overhang').selectedIndex;
    if (index === 13) {
        restr_update_rvcustom();
    } else {
        restr_update_rv();
    }
}


function restr_update_rvcustom() {
    var rv_overhang = $('#filler-rv-overhang').val() +
        $('#rvcustom-restr-site').val() +
        $('#scodon-rv-overhang').val();
    $('#rv-overhang').html(rv_overhang);
    mini_validator_rv(rv_overhang);
    end_notification_rv(rv_overhang);
}


function poly_update() {
    var polypeptides_labels = $('#polypeptides-labels');
    var polypeptides_viewer = $('#polypeptides-viewer');
    var i;

    polypeptides_labels.html('');
    polypeptides_viewer.html('');
    build_polypeptides();

    for (i = 0; i < CCD.polypeptides.length; i += 1) {
        polypeptides_labels.append('<label class="label-poly-check"><input type="checkbox" class="poly-check1">' +
            CCD.polypeptides[i][0] + '</label><br>');
        polypeptides_viewer.append('<span style="white-space:pre;">' + CCD.polypeptides[i][1] + '</span><br>');
    }
}


function select_all_serv() {
    $('.service-label > input').each(function () {
        this.checked = true;
    });
}


function deselect_all_serv() {
    $('.service-label > input').each(function () {
        this.checked = false;
    });
}


function plasmid_populate() {
    var polypeptides_list = $('#polypeptides-list');
    var i;

    build_polypeptides();
    polypeptides_list.html('');

    for (i = 0; i < CCD.polypeptides.length; i += 1) {
        polypeptides_list.append('<label class="label-poly-check"><input type="checkbox" class="poly-check-pm" id="">' +
            CCD.polypeptides[i][0] + '</label><br>');
    }
    $('#poly-check-select-pm')[0].checked = true;
    $('.poly-check-pm').each(function (i) {
        this.checked = true;
    });
    $.ajax({
        async: true,
        url: CCD.siteurl + 'get_separated_vectors',
        type: 'GET',
        success: function (response) {
            var json_obj = $.parseJSON(response);
            var vnames = json_obj.vnames;

            if (json_obj.status === 'VECTORS_ERROR') {
                alert('Error retrieving vectors data.\nPlease reload the page.');
                return false;
            }
            for (i = 0; i < vnames.length; i += 1) {
                $('#separated-vectors-list-pm').append('<option>' + vnames[i] + '</option>');
                $('#separated-vectors-list-pp').append('<option>' + vnames[i] + '</option>');
            }
        },
        error: function () {
            alert('Communication problem with CCD 2 server.\nPlease reload the page.');
        }
    });
}


function selectallnone_pm() {
    if ($('#poly-check-select-pm')[0].checked) { // if already pressed -> hide, return
        $('.poly-check-pm').each(function () {
            this.checked = true;
        });
    } else {
        $('.poly-check-pm').each(function () {
            this.checked = false;
        });
    }
}


function selectallnone_pp() {
    if ($('#poly-check-select-pp')[0].checked) { // if already pressed -> hide, return
        $('.poly-check1').each(function () {
            this.checked = true;
        });
        $('.poly-check2').each(function () {
            this.checked = true;
        });
        $('.poly-check3').each(function () {
            this.checked = true;
        });
    } else {
        $('.poly-check1').each(function () {
            this.checked = false;
        });
        $('.poly-check2').each(function () {
            this.checked = false;
        });
        $('.poly-check3').each(function () {
            this.checked = false;
        });
    }
}


function save_plasmid_maps() {
    var pname = '';
    var pstarts = '';
    var pstops = '';
    var pmfname = $('#pm-filename').val();
    var vname = $('#separated-vectors-list-pm option:selected').text();

    if (pmfname === null || pmfname === '') {
        pmfname = 'ccd_plasmid_maps.zip';
    } else {
        pmfname += '.zip';
    }
    $('.poly-check-pm').each(function (i) {
        if (this.checked) {
            pstarts += '&pstarts[]=' + CCD.polypeptides[i][2];
            pstops += '&pstops[]=' + CCD.polypeptides[i][3];
        }
    });
    if (pstarts === null || pstarts === '') {
        alert('Please select polypeptides.');
        return false;
    }
    if ($('#pp-name').val() === '') {
        pname = 'Pp';
    } else {
        pname = $('#pp-name').val();
    }
    var get_link_uri = CCD.siteurl + 'mapmake?vname=' + vname + '&pname=' + pname +
        '&pmfname=' + pmfname + pstarts + pstops;
    var popout = window.open(get_link_uri, '_blank');
    window.setTimeout(function () {
        popout.close();
    }, 1000);
}


function get_crysol() {
    if (CCD.crysol_task_ids) {
        crysol_revoke();
        setTimeout(crysol_requests, 2500);
    } else {
        crysol_requests();
    }
}


function save_crysol() {
    var link = document.createElement('a');
    var mimeType = 'text/plain';
    var pfname = $('#crysol-filename').val();
    var formatted_predictions = '';
    var labels_lines, secret_lines, proso_lines, proso2_lines, rpsp_lines;

    if (pfname === null || pfname === '') {
        pfname = 'ccd_crysol.csv';
    } else {
        pfname += '.csv';
    }
    link.setAttribute('download', pfname);

    if ($('#crysol-pp').css("display") !== 'none') {
        labels_lines = $('#crysol-pp').html().split('<br>');
        labels_lines.pop();
        secret_lines = $('#pp-secret-viewer').html().split('<br>');
        proso_lines = $('#pp-proso-viewer').html().split('<br>');
        proso2_lines = $('#pp-proso2-viewer').html().split('<br>');
        rpsp_lines = $('#pp-rpsp-viewer').html().split('<br>');

        jQuery.each(labels_lines, function (i) {
            formatted_predictions += labels_lines[i] + ',';
            formatted_predictions += secret_lines[i] + ',';
            formatted_predictions += proso_lines[i] + ',';
            formatted_predictions += proso2_lines[i] + ',';
            formatted_predictions += rpsp_lines[i] + '\n';
        });
    }
    if ($('#tagged-crysol-title').css("display") !== 'none') {
        labels_lines = $('#crysol-tagged').html().split('<br>');
        labels_lines.pop();
        secret_lines = $('#tagged-secret-viewer').html().split('<br>');
        proso_lines = $('#tagged-proso-viewer').html().split('<br>');
        proso2_lines = $('#tagged-proso2-viewer').html().split('<br>');
        rpsp_lines = $('#tagged-rpsp-viewer').html().split('<br>');

        jQuery.each(labels_lines, function (i) {
            formatted_predictions += labels_lines[i] + '(tagged),';
            formatted_predictions += secret_lines[i] + ',';
            formatted_predictions += proso_lines[i] + ',';
            formatted_predictions += proso2_lines[i] + ',';
            formatted_predictions += rpsp_lines[i] + '\n';
        });
    }
    if ($('#cleaved-crysol-title').css("display") !== 'none') {
        labels_lines = $('#crysol-cleaved').html().split('<br>');
        labels_lines.pop();
        secret_lines = $('#cleaved-secret-viewer').html().split('<br>');
        proso_lines = $('#cleaved-proso-viewer').html().split('<br>');
        proso2_lines = $('#cleaved-proso2-viewer').html().split('<br>');
        rpsp_lines = $('#cleaved-rpsp-viewer').html().split('<br>');

        jQuery.each(labels_lines, function (i) {
            formatted_predictions += labels_lines[i] + '(cleaved),';
            formatted_predictions += secret_lines[i] + ',';
            formatted_predictions += proso_lines[i] + ',';
            formatted_predictions += proso2_lines[i] + ',';
            formatted_predictions += rpsp_lines[i] + '\n';
        });
    }
    link.setAttribute('href', 'data:' + mimeType + ';charset=utf-8,' + encodeURIComponent(formatted_predictions));
    document.body.appendChild(link); // Firefox requires the link to be in the body
    link.click();
    document.body.removeChild(link); // remove the link when done
}


function realign() {
    raw_predictions();

    var flags_dictionary = {};

    if ($('#model-organisms')[0].checked === true) {
        flags_dictionary.model = true;
    } else {
        flags_dictionary.model = false;
    }

    if ($('#paralog-species')[0].checked === true) {
        flags_dictionary.within_species_paralog = true;
    } else {
        flags_dictionary.within_species_paralog = false;
    }

    if ($('#ortho-one2many')[0].checked === true) {
        flags_dictionary.ortholog_one2many = true;
    } else {
        flags_dictionary.ortholog_one2many = false;
    }

    if ($('#ortho-many2many')[0].checked === true) {
        flags_dictionary.ortholog_many2many = true;
    } else {
        flags_dictionary.ortholog_many2many = false;
    }

    if ($('#ortho-one2one')[0].checked === true) {
        flags_dictionary.ortholog_one2one = true;
    } else {
        flags_dictionary.ortholog_one2one = false;
    }

    $('#evolutionator-labels').html('');
    $('#alignments-alerts').html('');
    $('#evolutionator-viewer').html('Waiting for multiple sequence alignment...');
    $('#evolutionator-viewer').css('overflow-x', 'auto');
    $('#alignment-options').css('display', 'none');

    $.ajax({
        async: true,
        type: 'POST',
        url: CCD.siteurl + 'realign',
        data: {
            flags_dictionary: JSON.stringify(flags_dictionary)
        },
        success: function (response) {
            var json_obj = $.parseJSON(response);
            var realignment = json_obj.realignment;
            var realignment_flags = json_obj.realignment_flags;

            $('#evolutionator-viewer').html(''); // clear wait message
            if (realignment && json_obj.status === 'OK') {
                $('#alignment-options').css('display', 'inline');
                if (realignment_flags.model === true) {
                    $('#model-organisms')[0].checked = true;
                } else {
                    $('#model-organisms')[0].checked = false;
                }
                if (realignment_flags.within_species_paralog === true) {
                    $('#paralog-species')[0].checked = true;
                } else {
                    $('#paralog-species')[0].checked = false;
                }
                if (realignment_flags.ortholog_one2many === true) {
                    $('#ortho-one2many')[0].checked = true;
                } else {
                    $('#ortho-one2many')[0].checked = false;
                }
                if (realignment_flags.ortholog_many2many === true) {
                    $('#ortho-many2many')[0].checked = true;
                } else {
                    $('#ortho-many2many')[0].checked = false;
                }
                if (realignment_flags.ortholog_one2one === true) {
                    $('#ortho-one2one')[0].checked = true;
                } else {
                    $('#ortho-one2one')[0].checked = false;
                }
                CCD.multiple_alignment = realignment;
                draw_malignment();
                draw_predictions();
            } else {
                $('#alignments-alerts').html('<div class="alert alert-warning"><a class="close" data-dismiss="alert">×</a><span>Warning: ' +
                    json_obj.reason + '</span></div>');
            }
        },
        error: function (jqXHR, exception) {
            if (exception === 'abort') {
                $('#alignments-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                    '<span>Alignment request aborted.</span></div>');
            } else {
                $('#alignments-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                    '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
            }
        }
    });
}


function clear_sequence() {
    document.getElementById('dna-orf').value = '';
    document.getElementById('file-seq').type = '';
    document.getElementById('file-seq').type = 'file';
    $('#input-alerts').html('');
    $('.isoform-check').each(function (i) {
        this.checked = false;
    });
}


function save_alignment() {
    //preparing the alignment in fasta format
    var fasta = []
    var filename = ''
    var uniprot_id = $('#search-id').val();
    if (uniprot_id){
        filename = `${uniprot_id}_alignment.fasta`
    } else {
        filename = 'CCD_alignments.fasta'
    }
    for (var i = 0; i < CCD.multiple_alignment.length; i += 1){
        fasta.push(`>${CCD.multiple_alignment[i][0]}\n${CCD.multiple_alignment[i][1]}\n`);
    }
    var file = new Blob(fasta, {type: 'text/plain'});
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

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
    $('#alignments-alerts').html('');
    $('#evolutionator-labels').html('');
    $('#evolutionator-viewer').html('');
    $('#predictions-alerts-placeholder').html('');
    $('#predictions-info-placeholder').html('');
    $('sequences-labels').html('');
    $('sequences-viewer').html('');
    $('aa-label').html('');
    $('aa-viewer').html('');

    clear_after_primers();
}


function clear_after_primers() {
    $('#methods-overhangs-alerts-placeholder').html('');

    $('#primers-alerts-placeholder').html('');
    $('#primers-viewer').html('');

    // this is for plasmid maps panel
    $('#polypeptides-list').html('');

    $('#polypeptides-labels').html('');
    $('#polypeptides-viewer').html('');
    $('#MWs-viewer').html('');
    $('#pIs-viewer').html('');
    $('#epsilons-viewer').html('');

    clear_crysol();
}


function clear_crysol() {
    $('#polypeptides-alerts').html(''); // we need this alert, for no selected polypeptides to submit notification
    var i;

    if (CCD.crysol_task_ids) {
        crysol_revoke();
    }
    $('#crysol-info').html('');
    $('#crysol-errors').html('');

    clear_crysol_raw();
    clear_crysol_tagged();
    clear_crysol_cleaved();
}


function clear_crysol_raw() {
    $('#crysol-pp').html('');
    $('#pp-secret-viewer').html('');
    $('#pp-proso-viewer').html('');
    $('#pp-proso2-viewer').html('');
    $('#pp-rpsp-viewer').html('');
}


function clear_crysol_tagged() {
    $('#tagged-crysol-title').html('');
    $('#crysol-tagged').html('');
    $('#tagged-secret-viewer').html('');
    $('#tagged-proso-viewer').html('');
    $('#tagged-proso2-viewer').html('');
    $('#tagged-rpsp-viewer').html('');
}


function clear_crysol_cleaved() {
    $('#cleaved-crysol-title').html('');
    $('#crysol-cleaved').html('');
    $('#cleaved-secret-viewer').html('');
    $('#cleaved-proso-viewer').html('');
    $('#cleaved-proso2-viewer').html('');
    $('#cleaved-rpsp-viewer').html('');
}


function clear4tag() {
    $('#polypeptides-alerts').html('');
    if (CCD.tag_ajax) { // first time the variable should have a falshy value because it is declared as a global with a null value
        CCD.tag_ajax.abort(); // abort would fail in null variable
    }
    $('#tag-panel').html('');

}


function clear4cleave() {
    $('#polypeptides-alerts').html('');
    if (CCD.cleave_ajax) {
        CCD.cleave_ajax.abort();
    }
    $('#cleave-panel').html('');
}


// hide time
function hide_byclass(classname) {
    $('.' + classname).each(function () {
        $(this).css('display', 'none');
    });
}


// show time!
function show_byclass(classname) {
    $('.' + classname).each(function () {
        $(this).css('display', 'block');
    });
}


function poly_tag() {
    if ($('#poly-tagged-toggle').attr('aria-pressed') === 'true') { // if already pressed -> hide, return
        unpress_tag_button();
        clear4tag();
        hide_byclass('tag');
        return false;
    }

    show_byclass('tag');
    update_tagged();
    // Attention! These button-related operations can go below because the above ajax call in the function update_cleaved(vname, pname) is async: false = synchronous  
    $('#poly-tagged-toggle').html('Hide tagged constructs');
    $('#poly-tagged-toggle').addClass("active");
    $('#poly-tagged-toggle').attr("aria-pressed", "true");
}


function unpress_tag_button() {
    $('#poly-tagged-toggle').html('Show tagged constructs');
    $('#poly-tagged-toggle').removeClass("active");
    $('#poly-tagged-toggle').attr("aria-pressed", "false");
}


function poly_cleave() {
    if ($('#poly-cleaved-toggle').attr('aria-pressed') === 'true') { // if already pressed -> hide, return
        unpress_cleave_button();
        clear4cleave();
        hide_byclass('cleave');
        return false;
    }

    show_byclass('cleave');
    update_cleaved();
    // Attention! These button-related operations can go below because the above ajax call in the function update_cleaved(vname, pname) is async: false = synchronous  
    $('#poly-cleaved-toggle').html('Hide cleaved constructs');
    $('#poly-cleaved-toggle').addClass("active");
    $('#poly-cleaved-toggle').attr("aria-pressed", "true");
}


function unpress_cleave_button() {
    $('#poly-cleaved-toggle').html('Show cleaved constructs');
    $('#poly-cleaved-toggle').removeClass("active");
    $('#poly-cleaved-toggle').attr("aria-pressed", "false");
}


function sync_vector_pm() {
    var current_value = $('#separated-vectors-list-pm').val();
    $('#separated-vectors-list-pp').val(current_value);
    update_tagged();
    update_cleaved();
}


function sync_vector_pp() {
    var current_value = $('#separated-vectors-list-pp').val();
    $('#separated-vectors-list-pm').val(current_value);
    update_tagged();
    update_cleaved();
}


function update_tagged() {
    var vname = $('#separated-vectors-list-pp option:selected').text();
    var pname = '';

    if ($('#pp-name').val() === '') {
        pname = 'Pp';
    } else {
        pname = $('#pp-name').val();
    }

    clear4tag();
    $('#tag-info').html('Processing tagging...');

    CCD.tag_ajax = $.ajax({
        async: false,
        url: CCD.siteurl + 'tag_cleave',
        data: {
            polypeptides: JSON.stringify(CCD.polypeptides),
            vname: vname,
            pname: pname
        },
        type: 'POST',
        success: function (response) {
            var json_obj = $.parseJSON(response);
            CCD.tagged_cleaved = json_obj.tagged_cleaved;
            var tagged_pp_labels = '';
            var tagged_MWs = '';
            var tagged_pIs = '';
            var tagged_epsilons = '';
            var tagged_pp = '';
            var j;

            for (j = 0; j < CCD.tagged_cleaved.length; j += 1) {
                tagged_pp_labels += '<label class="label-poly-check"><input type="checkbox" class="poly-check2">' +
                    CCD.tagged_cleaved[j][0] + '</label><br>';
                tagged_MWs += CCD.tagged_cleaved[j][1] + '<br>';
                tagged_pIs += CCD.tagged_cleaved[j][2] + '<br>';
                tagged_epsilons += CCD.tagged_cleaved[j][3] + '<br>';
                tagged_pp += CCD.tagged_cleaved[j][4] + '<br>';
            }
            $('#tag-info').html('');
            $('#tag-panel').append('<br class="tag">' +
                '<br class="tag">' +
                '<h6 class="card-title tag" id="tagged-pp-title">Tagged constructs in vector ' + $('#separated-vectors-list-pp option:selected').text() + '</h6>' +
                '<table class="poly-panel tag">' +
                '<tr>' +
                '<th class="poly-labels"></th>' +
                '<th class="polysm-3 bold">MW</th>' +
                '<th class="polysm-3 bold">pI</th>' +
                '<th class="polysm-3 bold">ε</th>' +
                '</tr>' +
                '<tr>' +
                '<td class="poly-labels" id="">' + tagged_pp_labels + '</td>' +
                '<td class="poly-3" id="">' + tagged_MWs + '</td>' +
                '<td class="poly-3" id="">' + tagged_pIs + '</td>' +
                '<td class="poly-3" id="">' + tagged_epsilons + '</td>' +
                '<td><div class="poly-viewer" id="">' + tagged_pp + '</div></td>' +
                '</tr>' +
                '</table>');
        },
        error: function () {
            $('#polypeptides-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
        }
    });
}


function update_cleaved() {
    var vname = $('#separated-vectors-list-pp option:selected').text();
    var pname = '';

    if ($('#pp-name').val() === '') {
        pname = 'Pp';
    } else {
        pname = $('#pp-name').val();
    }

    clear4cleave();
    $('#cleave-info').html('Processing cleavage...');

    CCD.cleave_ajax = $.ajax({
        async: false,
        url: CCD.siteurl + 'tag_cleave',
        data: {
            polypeptides: JSON.stringify(CCD.polypeptides),
            vname: vname,
            pname: pname
        },
        type: 'POST',
        success: function (response) {
            var json_obj = $.parseJSON(response);
            CCD.tagged_cleaved = json_obj.tagged_cleaved;
            CCD.proteases = CCD.tagged_cleaved[0][9];
            var cleaved_pp_labels;
            var cleaved_MWs;
            var cleaved_pIs;
            var cleaved_epsilons;
            var cleaved_pp;
            var i, j;

            for (i = 0; i < CCD.proteases.length; i += 1) {
                cleaved_pp_labels = '';
                cleaved_MWs = '';
                cleaved_pIs = '';
                cleaved_epsilons = '';
                cleaved_pp = '';

                for (j = 0; j < CCD.tagged_cleaved.length; j += 1) {
                    cleaved_pp_labels += '<label class="label-poly-check"><input type="checkbox" class="poly-check3" name="' +
                        CCD.proteases[i] + '">' + CCD.tagged_cleaved[j][0] + '</label><br>';
                    cleaved_MWs += CCD.tagged_cleaved[j][5][i] + '<br>';
                    cleaved_pIs += CCD.tagged_cleaved[j][6][i] + '<br>';
                    cleaved_epsilons += CCD.tagged_cleaved[j][7][i] + '<br>';
                    cleaved_pp += '<span class="' + CCD.proteases[i] + '-cleaved">' + CCD.tagged_cleaved[j][8][i] + '</span><br>';
                }
                $('#cleave-info').html('');
                $('#cleave-panel').append('<br class="cleave">' +
                    '<h6 class="card-title cleave" id="cleaved-pp-title">Cleaved constructs with ' + CCD.proteases[i] + '</h6>' +
                    '<table class="poly-panel cleave">' +
                    '<tr>' +
                    '<th class="poly-labels"></th>' +
                    '<th class="polysm-3 bold">MW</th>' +
                    '<th class="polysm-3 bold">pI</th>' +
                    '<th class="polysm-3 bold">ε</th>' +
                    '</tr>' +
                    '<tr>' +
                    '<td class="poly-labels" id="">' + cleaved_pp_labels + '</td>' +
                    '<td class="poly-3" id="">' + cleaved_MWs + '</td>' +
                    '<td class="poly-3" id="">' + cleaved_pIs + '</td>' +
                    '<td class="poly-3" id="">' + cleaved_epsilons + '</td>' +
                    '<td><div class="poly-viewer" id="">' + cleaved_pp + '</div></td>' +
                    '</tr>' +
                    '</table>');
            }
        },
        error: function () {
            $('#polypeptides-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a><span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
        }
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


function draw_malignment() {
    var gaps_counter = 0;
    var gaps_start = 0;
    var gaps_spot = [];
    var position = 0;
    var gaps_insert = '';
    var new_prediction = '';
    var i, j, k, n;
    var painte_query = []; // label and painted sequence of query seq
    var multiple_alignment = CCD.multiple_alignment;
    var painte_malignment = paint_malignment(multiple_alignment);

    $('#evolutionator-labels').html('');
    $('#evolutionator-viewer').html('');
    for (i = 0; i < multiple_alignment.length; i += 1) {
        if (multiple_alignment[i][0].slice(-5) === 'Query') {
            CCD.query_str = multiple_alignment[i][1].replace(/-/g, ' ');
            painte_query[0] = multiple_alignment[i][0];
            painte_query[1] = painte_malignment[i];
        } else {
            $('#evolutionator-labels').append('<span class="malignment_seq' + i + '"><span onclick="minus_one(this)" style="cursor: pointer;">&#10134; </span>' +
                multiple_alignment[i][0] + '</span><br class="malignment_seq' + i + '">');
            $('#evolutionator-viewer').append('<span style="white-space:pre;" class="malignment_seq' + i + '">' +
                painte_malignment[i] + '</span><br class="malignment_seq' + i + '">');
        }
    }
    $('#evolutionator-viewer').css('overflow-x', 'scroll');
    // append the query in the end and indexes etc.
    $('#evolutionator-labels').append('<span class="malignment_query"><span onclick="minus_one(this)" style="cursor: pointer;">&#10134; </span>' +
        painte_query[0] + '</span><br class="malignment_query">');
    $('#evolutionator-viewer').append('<span style="white-space:pre;" class="malignment_query">' +
        painte_query[1] + '</span><br class="malignment_query">'); // append at the end
}


function paint_malignment(malignment) { // also replaces - with gap
    var row = '';
    var row_span = '';
    var column = '';
    var residue_group1_percentage = 0;
    var residue_group2_percentage = 0;
    var residue_group3_percentage = 0;
    var residue_group4_percentage = 0;
    var residue_group5_percentage = 0;
    var painte_malignment = [];
    var i, j, k;

    for (i = 0; i < malignment.length; i += 1) {
        row = malignment[i][1].replace(/-/g, ' ');
        row_span = '';
        for (j = 0; j < row.length; j += 1) {
            // core of painting

            if (row.charAt(j) === 'P') {
                row_span += '<span style="cursor:default; background-color: rgb(192, 192, 0);">' + row.charAt(j) + '</span>';
            } else if (row.charAt(j) === 'G') {
                row_span += '<span style="cursor:default; background-color: rgb(240, 144, 72);">' + row.charAt(j) + '</span>';
            } else if (row.charAt(j) === 'C') {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }
                residue_group1_percentage = (column.match(/[CS]/g) || []).length / column.length;
                residue_group2_percentage = (column.match(/[WLVIMAFCHP]/g) || []).length / column.length;

                if (residue_group1_percentage > 0.85) {
                    row_span += '<span style="cursor:default; background-color: rgb(128, 160, 240);">' + row.charAt(j) + '</span>';
                } else if (residue_group2_percentage > 0.6) {
                    row_span += '<span style="cursor:default; background-color: rgb(240, 128, 128);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else if (/[AILMFWV]/.test(row.charAt(j))) {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }

                if (row.charAt(j) === 'A') { // two rules for A, watch out!
                    residue_group1_percentage = (column.match(/[TSG]/g) || []).length / column.length;

                    if (residue_group1_percentage > 0.75) {
                        row_span += '<span style="cursor:default; background-color: rgb(128, 160, 240);">' + row.charAt(j) + '</span>';
                        continue;
                    }
                }

                residue_group1_percentage = (column.match(/[AWLVIMFCHPY]/g) || []).length / column.length;

                if (residue_group1_percentage > 0.6) {
                    row_span += '<span style="cursor:default; background-color: rgb(128, 160, 240);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else if (/[KR]/.test(row.charAt(j))) {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }
                residue_group1_percentage = (column.match(/[KR]/g) || []).length / column.length;
                residue_group2_percentage = (column.match(/[KRQ]/g) || []).length / column.length;
                residue_group3_percentage = (column.match(/[QN]/g) || []).length / column.length;

                if (residue_group1_percentage > 0.6 || residue_group2_percentage > 0.8 || residue_group3_percentage > 0.6) {
                    row_span += '<span style="cursor:default; background-color: rgb(240, 21, 5);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else if (/[E]/.test(row.charAt(j))) {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }
                residue_group1_percentage = (column.match(/[KR]/g) || []).length / column.length;
                residue_group2_percentage = (column.match(/[QE]/g) || []).length / column.length;
                residue_group3_percentage = (column.match(/[EQD]/g) || []).length / column.length;

                if (residue_group1_percentage > 0.6 || residue_group2_percentage > 0.5 || residue_group3_percentage > 0.85) {
                    row_span += '<span style="cursor:default; background-color: rgb(192, 72, 192);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else if (/[D]/.test(row.charAt(j))) {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }
                residue_group1_percentage = (column.match(/[KR]/g) || []).length / column.length;
                residue_group2_percentage = (column.match(/[KRQ]/g) || []).length / column.length;
                residue_group3_percentage = (column.match(/[ED]/g) || []).length / column.length;
                residue_group4_percentage = (column.match(/[N]/g) || []).length / column.length;

                if (residue_group1_percentage > 0.6 || residue_group2_percentage > 0.85 || residue_group3_percentage > 0.50 || residue_group4_percentage > 0.50) {
                    row_span += '<span style="cursor:default; background-color: rgb(192, 72, 192);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else if (/[N]/.test(row.charAt(j))) {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }
                residue_group1_percentage = (column.match(/[N]/g) || []).length / column.length;
                residue_group2_percentage = (column.match(/[NY]/g) || []).length / column.length;
                residue_group3_percentage = (column.match(/[DE]/g) || []).length / column.length;
                residue_group4_percentage = (column.match(/[QTS]/g) || []).length / column.length;

                if (residue_group1_percentage > 0.5 || residue_group2_percentage > 0.85 || residue_group3_percentage > 0.85 || residue_group4_percentage > 0.85) {
                    row_span += '<span style="cursor:default; background-color: rgb(21, 192, 21);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else if (/[Q]/.test(row.charAt(j))) {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }
                residue_group1_percentage = (column.match(/[KR]/g) || []).length / column.length;
                residue_group2_percentage = (column.match(/[QE]/g) || []).length / column.length;
                residue_group3_percentage = (column.match(/[QEKR]/g) || []).length / column.length;
                residue_group4_percentage = (column.match(/[QN]/g) || []).length / column.length;
                residue_group5_percentage = (column.match(/[DE]/g) || []).length / column.length;

                if (residue_group1_percentage > 0.6 || residue_group2_percentage > 0.5 || residue_group3_percentage > 0.85 || residue_group4_percentage > 0.5 || residue_group5_percentage > 0.6) {
                    row_span += '<span style="cursor:default; background-color: rgb(21, 192, 21);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else if (/[ST]/.test(row.charAt(j))) {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }
                residue_group1_percentage = (column.match(/[WLVIMAFCHP]/g) || []).length / column.length;
                residue_group2_percentage = (column.match(/[TS]/g) || []).length / column.length;
                residue_group3_percentage = (column.match(/[QST]/g) || []).length / column.length; // unsure about this rule since the table was missing something there

                if (residue_group1_percentage > 0.6 || residue_group2_percentage > 0.5 || residue_group3_percentage > 0.85) {
                    row_span += '<span style="cursor:default; background-color: rgb(21, 192, 21);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else if (/[HY]/.test(row.charAt(j))) {
                column = '';
                for (k = 0; k < malignment.length; k += 1) {
                    column += malignment[k][1][j];
                }

                if (row.charAt(j) === 'H') {
                    residue_group1_percentage = (column.match(/[KR]/g) || []).length / column.length;

                    if (residue_group1_percentage > 0.75) {
                        row_span += '<span style="cursor:default; background-color: rgb(21, 164, 164);">' + row.charAt(j) + '</span>';
                        continue;
                    }
                }

                residue_group1_percentage = (column.match(/[WLVIMAFCHP]/g) || []).length / column.length;
                residue_group2_percentage = (column.match(/[WYACPQFHILMV]/g) || []).length / column.length;

                if (residue_group1_percentage > 0.6 || residue_group2_percentage > 0.85) {
                    row_span += '<span style="cursor:default; background-color: rgb(21, 164, 164);">' + row.charAt(j) + '</span>';
                } else {
                    row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
                }
            } else {
                row_span += '<span style="cursor:default;">' + row.charAt(j) + '</span>';
            }
        }
        painte_malignment[i] = row_span;
    }
    return painte_malignment;
}


function draw_predictions() {
    var alignment_visibility = $('#alignment-toggler').is(':checked');
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
    $('.service-label').each(function () {
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
                                                            CCD.gapped_predictions[n][1].slice(position)].join('');
                            if (CCD.gapped_predictions[n][2]) { // some predictions have an extra element, which also needs to be gapped
                                CCD.gapped_predictions[n][2] = [CCD.gapped_predictions[n][2].slice(0, position), 
                                                                gaps_insert, CCD.gapped_predictions[n][2].slice(position)].join('');
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
                    for (var n = 0; n < CCD.gapped_predictions.length; n += 1){ 
                        if (CCD.gapped_predictions[n][0] === 'PDB_95'){
                            //first we need to map gaps_spot back to the sequence's number. 
                            var seq = CCD.gapped_predictions[n][1];
                            var structures = JSON.parse(CCD.gapped_predictions[n][3]);
                            var positions = Object.keys(structures)
                            var counter = 0;
                            for (var i = 0; i <= seq.length; i += 1){
                                if (seq[i] === '>' || seq[i] === '<' || seq[i] === String.fromCharCode('9674')){
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
                        case 'NLS':
                            draw_nls(i);
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
            }
            else {
                document.getElementById("sequences-labels").style.paddingBottom = "0px";
            }
        }
    });
}


function build_polypeptides() {
    var vstarts = CCD.vstarts;
    var vstops = CCD.vstops;
    CCD.polypeptides = [];
    var polypeptide_sequence = '';
    var polypeptide_name = '';
    var gaps_reduced;
    var gaps;
    var raw_poly;
    var i, j;

    for (i = 0; i < vstarts.length; i += 1) {
        for (j = 0; j < vstops.length; j += 1) {
            if (vstarts[i] < vstops[j]) {
                if ($('#pp-name').val() === '') {
                    polypeptide_name = 'Pp_' + vstarts[i] + '_' + vstops[j];
                } else {
                    polypeptide_name = $('#pp-name').val() + '_' + vstarts[i] + '_' + vstops[j];
                }
                if ($('#poly-align-toggle')[0].checked) {
                    gaps = '';
                    if (i !== 0) {
                        gaps_reduced = parseInt(vstarts[i]) - parseInt(vstarts[0]);
                        gaps = Array(gaps_reduced + 1).join(' ');
                    }
                    raw_poly = CCD.AAseq.slice(vstarts[i] - 1, vstops[j]);
                    polypeptide_sequence = gaps.concat(raw_poly);
                } else {
                    polypeptide_sequence = CCD.AAseq.slice(vstarts[i] - 1, vstops[j]);
                }
            } else {
                continue;
            }
            CCD.polypeptides.push([polypeptide_name, polypeptide_sequence, vstarts[i], vstops[j]]);
        }
    }
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


function mini_validator_fw(fw) {
    // clear previous warning messages
    $('#methods-fw-overhang-warnings-placeholder').html('');
    //Search for charaters that are not G, A, T or C.
    if (fw.search(/^[AaTtGgCc]*$/g) === -1) {
        //The seq string contains non-DNA characters
        $('#methods-fw-overhang-warnings-placeholder').html('<div class="alert alert-warning"><a class="close" data-dismiss="alert">×</a>' +
            '<span>Detected non DNA characters in forward overhang.</span></div>');
    }
}


function mini_validator_rv(rv) {
    // clear previous warning messages
    $('#methods-rv-overhang-warnings-placeholder').html('');
    //Search for charaters that are not G, A, T or C.
    if (rv.search(/^[AaTtGgCc]*$/g) === -1) {
        //The seq string contains non-DNA characters
        $('#methods-rv-overhang-warnings-placeholder').html('<div class="alert alert-warning"><a class="close" data-dismiss="alert">×</a>' +
            '<span>Detected non DNA characters in reverse overhang.</span></div>');
    }
}


function end_notification_fw(fw) {
    // clear previous warning messages
    $('#fw-overhang-alerts').html('');
    if (fw.slice(-3) === 'atg') {
        $('#fw-overhang-alerts').html('<span style="display: inline-block; font-weight: bold;">Forward overhang ends with a start codon (atg).</span>');
        $('#fw-overhang').html(fw.slice(0, -3) + '<span style="color: red;">atg</span>');
    }
}


function end_notification_rv(rv) {
    // clear previous warning messages
    $('#rv-overhang-alerts').html('');
    if (rv.slice(-3) === 'tta' || rv.slice(-3) === 'tca' || rv.slice(-3) === 'cta') {
        $('#rv-overhang-alerts').html('<span style="display: inline-block; font-weight: bold;">Reverse overhang ends with a stop codon (' + rv.slice(-3) + ').</span>');
        $('#rv-overhang').html(rv.slice(0, -3) + '<span style="color: red;">' + rv.slice(-3) + '</span>');
    }
}


function match_vectors() {
    var bac_vectors_sublist = '';
    var ins_vectors_sublist = '';
    var mam_vectors_sublist = '';
    var i;
    var vectors_sublist = '';
    var max_len = Math.max(CCD.bac_vectors.length, CCD.ins_vectors.length, CCD.mam_vectors.length);
    var fw_overhang = $('#fw-overhang').text();
    var rv_overhang = $('#rv-overhang').text();

    for (i = 0; i < max_len; i += 1) {
        if (i < CCD.bac_vectors.length && fw_overhang === CCD.bac_vectors[i].fw_overhang && rv_overhang === CCD.bac_vectors[i].rv_overhang) {
            bac_vectors_sublist += CCD.bac_vectors[i].name + '<br>';
        }
        if (i < CCD.ins_vectors.length && fw_overhang === CCD.ins_vectors[i].fw_overhang && rv_overhang === CCD.ins_vectors[i].rv_overhang) {
            ins_vectors_sublist += CCD.ins_vectors[i].name + '<br>';
        }
        if (i < CCD.mam_vectors.length && fw_overhang === CCD.mam_vectors[i].fw_overhang && rv_overhang === CCD.mam_vectors[i].rv_overhang) {
            mam_vectors_sublist += CCD.mam_vectors[i].name + '<br>';
        }
    }
    if (bac_vectors_sublist !== '') {
        vectors_sublist += 'Bacterial<br>' + bac_vectors_sublist;
    }
    if (ins_vectors_sublist !== '') {
        vectors_sublist += 'Insect<br>' + ins_vectors_sublist;
    }
    if (mam_vectors_sublist !== '') {
        vectors_sublist += 'Mammalian<br>' + mam_vectors_sublist;
    }
    $('#matching-vectors').html(vectors_sublist);
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


function draw_pdb(ind) {
    var pred = CCD.gapped_predictions[ind][0];
    var symbols = CCD.gapped_predictions[ind][1];
    var hit_freq = CCD.gapped_predictions[ind][2];
    var structures = JSON.parse(CCD.gapped_predictions[ind][3]);
    //adding tooltips to the prediction label
    if (pred === 'PDB_95') {
        var tooltip_text = `Indicates region of the protein for which a >95% identical structure is known. 
    Darker shades of blue indicate progressively large number of structures that span the region.
    '>' indicates that a construct corresponding to a PDB entry starts around that specific residue*;
    '<' indicates that a construct corresponding to PDB entry stops at that residue*.
    '\u25CA' indicates a simultaneous start and stop.\n
    Hover on the start / stop to see which PDB structures match. Format: <PDBID>_<CHAIN>\n
    *: due to the presence of tags, BLAST might be confused about the exact start; you can always check the PDB id and verify the exact boundaries online.`;
        var spacer = "&nbsp".repeat(3)
    } else if (pred === 'PDB_50to95') {
        var tooltip_text = `BLAST against PDB - near homology.
Color indicates the region spanned by blast high scoring pairs with similarity between 100% and 50%.
Darker shades of blue indicate a larger number of individual hits in the database.`;
        var spacer = "&nbsp".repeat(3)
    } else if (pred === 'PDB_30to50') {
        var tooltip_text = `BLAST against PDB - remote homology
Color indicates the region spanned by blast high scoring pairs with similarity between 100% and 30%.
Darker shades of blue indicate a larger number of individual hits in the database`;
        var spacer = "&nbsp".repeat(3)
    }
    var tooltip_html = `<span data-toggle="tooltip" title="${tooltip_text}" name="2">${spacer}[?]</span>`
    $('#sequences-labels').append('<span>' + pred + tooltip_html + '</span><br>');

    // decorate by color and add starts stops if needed for PDB_95
    var pdb_span = '';
    var aa_tooltip_text = '';
    var flag = 1;
    var i, j, k, s_list, bg_color;
    var fg_color = 'rgb(0, 0, 0)';
    var bg_color = 'rgb(255, 255, 255)';
    for (i = 0; i < hit_freq.length; i += 1) {

        // for starts and stops of PDB_95 prediction, we want to add a popup that lists PDB structure that have that start/stop
        // data comes in as JSON of format {"aa_position": [STR1_CHAIN, STR2_CHAIN, ..]}
        if (pred === 'PDB_95') {
            k = i.toString();
            s_list = ''
            if (structures.hasOwnProperty(k)) {
                for (j = 0; j < structures[k].length; j += 1) {
                    s_list = structures[k].join('\n');
                    aa_tooltip_text = `data-toggle="tooltip" title="${s_list}"`;
                }
            } else {
                aa_tooltip_text = '';
            }
        }
        // standard coloring stuff
        if (hit_freq[i] === '1') {
            bg_color = 'rgb(230, 230, 255)';
        } else if (hit_freq[i] === '2' || hit_freq[i] === '3') {
            bg_color = 'rgb(204, 204, 255)';
        } else if (hit_freq[i] === '4' || hit_freq[i] === '5') {
            bg_color = 'rgb(179, 179, 255)';
        } else if (hit_freq[i] === '6' || hit_freq[i] === '7') {
            bg_color = 'rgb(153, 153, 255)';
        } else if (hit_freq[i] === '8' || hit_freq[i] === '9') {
            bg_color = 'rgb(128, 128, 255)';
        } else if (hit_freq[i] === '0') {
            bg_color = 'rgb(43, 43, 255)';
            fg_color = 'rgb(255, 255, 255)';
        } else {
            bg_color = 'rgb(255, 255, 255)'
        }
        //put tooltip text last or we lose the background color for some reason
        pdb_span += `<span style="cursor:default; background-color:${bg_color}; color:${fg_color} ${aa_tooltip_text}">${symbols[i]}</span>`;
    }
    $('#sequences-viewer').append('<span style="white-space:pre;">' + pdb_span + '</span><br>');
}


function draw_nls(ind) {
    $('#sequences-labels').append('<span>' + CCD.gapped_predictions[ind][0] + '</span><br>');
    var nls_pred = CCD.gapped_predictions[ind][1];
    var nls_span = '';
    var i;

    for (i = 0; i < nls_pred.length; i += 1) {
        if (nls_pred[i] === '5') {
            nls_span += '<span style="cursor:default; background-color:rgb(230, 230, 255);">' + nls_pred[i] + '</span>';
        } else if (nls_pred[i] === '6') {
            nls_span += '<span style="cursor:default; background-color:rgb(204, 204, 255);">' + nls_pred[i] + '</span>';
        } else if (nls_pred[i] === '7') {
            nls_span += '<span style="cursor:default; background-color:rgb(179, 179, 255);">' + nls_pred[i] + '</span>';
        } else if (nls_pred[i] === '8') {
            nls_span += '<span style="cursor:default; background-color:rgb(153, 153, 255);">' + nls_pred[i] + '</span>';
        } else if (nls_pred[i] === '9') {
            nls_span += '<span style="cursor:default; background-color:rgb(128, 128, 255);">' + nls_pred[i] + '</span>';
        } else if (nls_pred[i] === '0') {
            nls_span += '<span style="cursor:default; background-color:rgb(43, 43, 255); color: rgb(255, 255, 255)">' + nls_pred[i] + '</span>';
        } else {
            nls_span += '<span style="cursor:default;">' + nls_pred[i] + '</span>';
        }
    }
    $('#sequences-viewer').append('<span style="white-space:pre;">' + nls_span + '</span><br>');
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


function draw_prediction(ind) {
    $('#sequences-labels').append('<span>' + CCD.predictions[ind][0] + '</span><br>');
    $('#sequences-viewer').append(CCD.predictions[ind][1] + '<br>');
}


function draw_wait(service_name) {
    $('#sequences-labels').append(service_name + '<br>');
    $('#sequences-viewer').append('Processing...<br>');
}


function check_pending() {
    var pending_flag = false; // preassume that every task has finished
    var i, j;

    for (i = 0; i < 3; i += 1) {
        for (j = 0; j < 3; j += 1) {
            if (CCD.crysol_pending[i][j] !== '' && CCD.crysol_pending[i][j] === true) {
                pending_flag = true;
            }
        }
    }
    if (pending_flag === false) {
        $('#crysol-info').html('');
    }
}


function rpsp_task(polypeptides_sequences, crysol_matrix, pp_afterfix) {
    $.ajax({
        type: 'POST',
        url: CCD.siteurl + 'crysol_prediction/RPSP',
        data: {
            selected_polypeptides_sequences: JSON.stringify(polypeptides_sequences)
        },
        async: true,
        success: function (data, status, request) {
            var status_url = request.getResponseHeader('Location');
            var task_id = status_url.substr(status_url.lastIndexOf('/') + 1);
            CCD.crysol_task_ids.push(task_id);
            update_rpsp(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix);
        },
        error: function () {
            $('#crysol-errors').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
            $('#crysol-info').html('');
        }
    });
}


function update_rpsp(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix) {
    if (CCD.crysol_stop_recursion) {
        return;
    }
    // send GET request to status URL
    $.getJSON(status_url, function (data) {
        if (data['state'] === 'SUCCESS') {
            var crysol_errors = $('#crysol-errors');
            var formatted_errors = '';
            var rpsp = data['result'];
            var i;
            var pending_index = '';

            if (pp_afterfix === ' (raw)') {
                pending_index = 0;
            } else if (pp_afterfix === ' (tagged)') {
                pending_index = 1;
            } else if (pp_afterfix === ' (cleaved)') {
                pending_index = 2;
            }
            CCD.crysol_pending[pending_index][0] = false;

            if (rpsp === 'ERROR' || rpsp === undefined) {
                formatted_errors += 'RPSP' + pp_afterfix + ' bad results or server unaccessible<br>';
                rpsp = new Array(polypeptides_sequences.length).fill('undefined');
            }
            for (i = 0; i < rpsp.length; i += 1) {
                crysol_matrix[i][1] = rpsp[i];
            }
            check_pending();
            sort_draw_crysol(pp_afterfix);
            if (formatted_errors !== '') {
                crysol_errors.append('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a><span>' +
                    formatted_errors + '</span></div>');
            }
        } else {
            // rerun in 2 seconds
            setTimeout(function () {
                update_rpsp(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix);
            }, 2000);
        }
    });
}


function proso2_task(polypeptides_sequences, crysol_matrix, pp_afterfix) {
    $.ajax({
        type: 'POST',
        url: CCD.siteurl + 'crysol_prediction/PROSO2',
        data: {
            selected_polypeptides_sequences: JSON.stringify(polypeptides_sequences)
        },
        async: true,
        success: function (data, status, request) {
            var status_url = request.getResponseHeader('Location');
            var task_id = status_url.substr(status_url.lastIndexOf('/') + 1);
            CCD.crysol_task_ids.push(task_id);
            update_proso2(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix);
        },
        error: function () {
            $('#crysol-errors').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
            $('#crysol-info').html('');
        }
    });
}


function update_proso2(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix) {
    if (CCD.crysol_stop_recursion) {
        return;
    }
    // send GET request to status URL
    $.getJSON(status_url, function (data) {
        if (data['state'] === 'SUCCESS') {
            var crysol_errors = $('#crysol-errors');
            var formatted_errors = '';
            var proso2 = data['result'];
            var i;
            var pending_index = '';

            if (pp_afterfix === ' (raw)') {
                pending_index = 0;
            } else if (pp_afterfix === ' (tagged)') {
                pending_index = 1;
            } else if (pp_afterfix === ' (cleaved)') {
                pending_index = 2;
            }
            CCD.crysol_pending[pending_index][1] = false;

            if (proso2 === 'ERROR' || proso2 === undefined) {
                formatted_errors += 'PROSOII' + pp_afterfix + ' bad results or server unaccessible<br>';
                proso2 = new Array(polypeptides_sequences.length).fill('undefined');
            } else {
                for (i = 0; i < polypeptides_sequences.length; i += 1) {
                    if (polypeptides_sequences[i].length < 21 || polypeptides_sequences[i].length > 2000) {
                        formatted_errors += crysol_matrix[i][0] + pp_afterfix + ', ';
                        proso2[i] = 'undefined';
                    }
                }
                if (formatted_errors !== '') {
                    formatted_errors = formatted_errors.substring(0, formatted_errors.length - 2);
                    formatted_errors = 'PROSOII only allows sequences of length between 21-2000, polypeptides with name ' + formatted_errors + ' do not comply.<br>'
                }
            }
            for (i = 0; i < proso2.length; i += 1) {
                crysol_matrix[i][2] = proso2[i];
            }
            check_pending();
            sort_draw_crysol(pp_afterfix);
            if (formatted_errors !== '') {
                crysol_errors.append('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a><span>' +
                    formatted_errors + '</span></div>');
            }
        } else {
            // rerun in 2 seconds
            setTimeout(function () {
                update_proso2(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix);
            }, 2000);
        }
    });
}


function secret_task(polypeptides_sequences, crysol_matrix, pp_afterfix) {
    $.ajax({
        type: 'POST',
        url: CCD.siteurl + 'crysol_prediction/SECRET',
        data: {
            selected_polypeptides_sequences: JSON.stringify(polypeptides_sequences)
        },
        async: true,
        success: function (data, status, request) {
            var status_url = request.getResponseHeader('Location');
            var task_id = status_url.substr(status_url.lastIndexOf('/') + 1);
            CCD.crysol_task_ids.push(task_id);
            update_secret(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix);
        },
        error: function () {
            $('#crysol-errors').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>' +
                '<span>Communication problem with CCD 2 server.<br>Please try again or reload the page.</span></div>');
            $('#crysol-info').html('');
        }
    });
}


function update_secret(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix) {
    if (CCD.crysol_stop_recursion) {
        return;
    }
    // send GET request to status URL
    $.getJSON(status_url, function (data) {
        if (data['state'] === 'SUCCESS') {
            var crysol_errors = $('#crysol-errors');
            var formatted_errors = '';
            var secret = data['result'];
            var i;
            var pending_index = '';
            if (pp_afterfix === ' (raw)') {
                pending_index = 0;
            } else if (pp_afterfix === ' (tagged)') {
                pending_index = 1;
            } else if (pp_afterfix === ' (cleaved)') {
                pending_index = 2;
            }
            CCD.crysol_pending[pending_index][2] = false;

            if (secret === 'ERROR' || secret === undefined) {
                formatted_errors += 'SECRET' + pp_afterfix + ' bad results or server unaccessible<br>';
                secret = new Array(polypeptides_sequences.length).fill('undefined');
            } else {
                for (i = 0; i < polypeptides_sequences.length; i += 1) {
                    if (polypeptides_sequences[i].length < 46 || polypeptides_sequences[i].length > 200) {
                        formatted_errors += crysol_matrix[i][0] + pp_afterfix + ', ';
                        secret[i] = 'undefined';
                    }
                }
                if (formatted_errors !== '') {
                    formatted_errors = formatted_errors.substring(0, formatted_errors.length - 2);
                    formatted_errors = 'SECRET only allows sequences of length between 46-200, polypeptides with name ' + formatted_errors + ' do not comply.<br>'
                }
            }
            for (i = 0; i < secret.length; i += 1) {
                crysol_matrix[i][3] = secret[i];
            }
            check_pending();
            sort_draw_crysol(pp_afterfix);
            if (formatted_errors !== '') {
                crysol_errors.append('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a><span>' +
                    formatted_errors + '</span></div>');
            }
        } else {
            // rerun in 2 seconds
            setTimeout(function () {
                update_secret(status_url, polypeptides_sequences, crysol_matrix, pp_afterfix);
            }, 2000);
        }
    });
}


function evolutionator_call() {
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

    // clear, hide
    $('#evolutionator-labels').html('');
    $('#alignments-alerts').html('');
    $('#evolutionator-viewer').html('Waiting for multiple sequence alignment...');
    $('#evolutionator-viewer').css('overflow-x', 'auto');
    $('.alignment-options').css('display', 'none');

    $.ajax({
        async: true,
        type: 'POST',
        url: CCD.siteurl + 'evolutionator',
        data: {
            uniprot_id: $('#search-id').val().trim(),
            AAseq: CCD.AAseq,
            flags_dictionary: JSON.stringify(flags_dictionary)
        },
        success: function (response) {
            var json_obj = $.parseJSON(response);
            var multiple_alignment = json_obj.multiple_alignment;
            var flags = json_obj.alignment_flags;
            var hits_by_flag = json_obj.hits_by_flag;

            $('#evolutionator-viewer').html(''); // clear wait message
            if (multiple_alignment && json_obj.status === 'OK') {
                $('.alignment-options').css('display', 'inline');
                if (flags.model === true) {
                    $('#model-organisms')[0].checked = true;
                } else {
                    $('#model-organisms')[0].checked = false;
                }
                $('#model-organisms').parent().append(' (' + hits_by_flag['model'] + ')');
                if (flags.within_species_paralog === true) {
                    $('#paralog-species')[0].checked = true;
                } else {
                    $('#paralog-species')[0].checked = false;
                }
                $('#paralog-species').parent().append(' (' + hits_by_flag['within_species_paralog'] + ')');
                if (flags.ortholog_one2many === true) {
                    $('#ortho-one2many')[0].checked = true;
                } else {
                    $('#ortho-one2many')[0].checked = false;
                }
                $('#ortho-one2many').parent().append(' (' + hits_by_flag['ortholog_one2many'] + ')');
                if (flags.ortholog_many2many === true) {
                    $('#ortho-many2many')[0].checked = true;
                } else {
                    $('#ortho-many2many')[0].checked = false;
                }
                $('#ortho-many2many').parent().append(' (' + hits_by_flag['ortholog_many2many'] + ')');
                if (flags.ortholog_one2one === true) {
                    $('#ortho-one2one')[0].checked = true;
                } else {
                    $('#ortho-one2one')[0].checked = false;
                }
                $('#ortho-one2one').parent().append(' (' + hits_by_flag['ortholog_one2one'] + ')');

                CCD.multiple_alignment = multiple_alignment;
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
    $('.service-label').each(function () {
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


function raw_malignment() {
    var query = [];
    var i;
    var multiple_alignment = CCD.multiple_alignment;

    $('#evolutionator-labels').html('');
    $('#evolutionator-viewer').html('');
    for (i = 0; i < multiple_alignment.length; i += 1) {
        if (multiple_alignment[i][0].slice(-5) === 'Query') {
            CCD.query_str = multiple_alignment[i][1].replace(/-/g, ' ');
            query[0] = multiple_alignment[i][0];
            query[1] = multiple_alignment[i][1];
        } else {
            $('#evolutionator-labels').append('<span class="malignment_seq' + i + '"><span onclick="minus_one(this)" style="cursor: pointer;">&#10134; </span>' +
                multiple_alignment[i][0] + '</span><br class="malignment_seq' + i + '">');
            $('#evolutionator-viewer').append('<span style="white-space:pre;" class="malignment_seq' + i + '">' +
                multiple_alignment[i][1].replace(/-/g, ' ') + '</span><br class="malignment_seq' + i + '">');
        }
    }
    $('#evolutionator-viewer').css('overflow-x', 'scroll');
    // append the query in the end and indexes etc.
    $('#evolutionator-labels').append('<span class="malignment_query"><span onclick="minus_one(this)" style="cursor: pointer;">&#10134; </span>' +
        query[0] + '</span><br class="malignment_query">');
    $('#evolutionator-viewer').append('<span style="white-space:pre;" class="malignment_query">' +
        query[1].replace(/-/g, ' ') + '</span><br class="malignment_query">'); // append at the end
}


function arrow_toggle(anchor) {
    var current_arrow = $(anchor).html();
    $(anchor).html(current_arrow === 'more' ? 'less' : 'more'); // php style if
}


function initialize_reset() {
    CCD.starts = [];
    CCD.stops = [];
    $('#ortho-one2one')[0].checked = true;
    $('#ortho-one2many')[0].checked = true;
    $('#ortho-many2many')[0].checked = false;
    $('#paralog-species')[0].checked = false;
    $('#model-organisms')[0].checked = true;

    $('#ortho-one2one').parent().contents().filter(function () {
        return this.nodeType === 3;
    }).remove();
    $('#ortho-one2one').parent().append(' 1 to 1 orthologs');
    $('#ortho-one2many').parent().contents().filter(function () {
        return this.nodeType === 3;
    }).remove();
    $('#ortho-one2many').parent().append(' 1 to many orthologs');
    $('#ortho-many2many').parent().contents().filter(function () {
        return this.nodeType === 3;
    }).remove();
    $('#ortho-many2many').parent().append(' many to many orthologs');
    $('#model-organisms').parent().contents().filter(function () {
        return this.nodeType === 3;
    }).remove();
    $('#model-organisms').parent().append(' representative species only');
    // paralogs to be removed
}


function draw_alignment_table(aligned_isoforms) {
    var formatted_labels = '';
    var formatted_aligned_seqs = '';
    var curr_length = 0;
    var indexes_string = '';
    var alignment_table = '';

    aligned_isoforms.sort(function (a, b) { // sort isoforms based on the name
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


function sort_draw_crysol(afterfix) {
    var i;

    if (afterfix === ' (raw)') {
        clear_crysol_raw();
        var sort_raw = CCD.crysol_raw.slice(); // deep copy the values of array

        show_byclass('pp-crysol');
        $('#crysol-save').css('display', 'block');
        $('#raw-crysol-title').html('Raw constructs');
        if (CCD.crysol_raw_sort_switch) {
            sort_constructs(sort_raw, CCD.crysol_raw_sort_switch);
        }
        for (i = 0; i < sort_raw.length; i += 1) {
            $('#crysol-pp').append(sort_raw[i][0] + '<br>');
            $('#pp-rpsp-viewer').append(sort_raw[i][1] + '<br>');
            $('#pp-proso2-viewer').append(sort_raw[i][2] + '<br>');
            $('#pp-secret-viewer').append(sort_raw[i][3] + '<br>');
        }
    }
    if (afterfix === ' (tagged)') {
        clear_crysol_tagged();
        var sort_tagged = CCD.crysol_tagged.slice(); // deep copy the values of the array

        show_byclass('tag-crysol');
        $('#crysol-save').css('display', 'block');
        $('#tagged-crysol-title').html('Tagged constructs in vector ' + $('#separated-vectors-list-pp option:selected').text());
        if (CCD.crysol_tagged_sort_switch) {
            sort_constructs(sort_tagged, CCD.crysol_tagged_sort_switch);
        }
        for (i = 0; i < sort_tagged.length; i += 1) {
            $('#crysol-tagged').append(sort_tagged[i][0] + '<br>');
            $('#tagged-rpsp-viewer').append(sort_tagged[i][1] + '<br>');
            $('#tagged-proso2-viewer').append(sort_tagged[i][2] + '<br>');
            $('#tagged-secret-viewer').append(sort_tagged[i][3] + '<br>');
        }
    }
    if (afterfix === ' (cleaved)') {
        clear_crysol_cleaved();
        var sort_cleaved = CCD.crysol_cleaved.slice(); // deep copy the values of the array

        show_byclass('cleave-crysol');
        $('#crysol-save').css('display', 'block');
        $('#cleaved-crysol-title').html('Cleaved constructs');
        if (CCD.crysol_cleaved_sort_switch) {
            sort_constructs(sort_cleaved, CCD.crysol_cleaved_sort_switch);
        }
        for (i = 0; i < sort_cleaved.length; i += 1) {
            $('#crysol-cleaved').append(sort_cleaved[i][0] + '<br>');
            $('#cleaved-rpsp-viewer').append(sort_cleaved[i][1] + '<br>');
            $('#cleaved-proso2-viewer').append(sort_cleaved[i][2] + '<br>');
            $('#cleaved-secret-viewer').append(sort_cleaved[i][3] + '<br>');
        }
    }
}


function sort_constructs(crysol_matrix, prediction_service) { // this operates by reference to the matrix
    var service2index = {
        'RPSP': 1,
        'PROSOII': 2,
        'SECRET': 3
    };

    crysol_matrix.sort(function (a, b) { // sort isoforms based on the name
        var keyA = a[service2index[prediction_service]],
            keyB = b[service2index[prediction_service]];
        // Compare the 2 keys
        if (keyA === 'undefined') {
            keyA = 0
        }
        if (keyB === 'undefined') {
            keyB = 0
        }

        if (parseFloat(keyA) > parseFloat(keyB)) return -1;
        if (parseFloat(keyA) < parseFloat(keyB)) return 1;
        return 0;
    });
}


function crysol_revoke() {
    CCD.crysol_stop_recursion = true; // all the recirsive functions check this variable before continuing..
    var task_ids_list = '';
    var i;

    for (i = 0; i < CCD.crysol_task_ids.length; i += 1) {
        task_ids_list += '&task_ids_list[]=' + CCD.crysol_task_ids[i];
    }
    var revokeall_tasks_uri = CCD.siteurl + 'revokeall?' + task_ids_list;

    $.ajax({
        type: 'GET',
        url: revokeall_tasks_uri,
        async: true,
        success: function () {},
        error: function () {
            alert('Could not revoke crystallisation and solubility tasks on server.');
        }
    });
}


function crysol_requests() {
    var selected_polypeptides_names = [];
    var tagged_selected_polypeptides_names = [];
    var cleaved_selected_polypeptides_names = [];
    var selected_polypeptides_sequences = [];
    var tagged_selected_polypeptides_sequences = [];
    var cleaved_selected_polypeptides_sequences = [];
    var vname = $('#separated-vectors-list-pp option:selected').text();
    var i, j;

    clear_crysol();
    CCD.crysol_pending = [['', '', ''], ['', '', ''], ['', '', '']]; // initialize pending list because is declared as null in globals array
    CCD.crysol_task_ids = [];
    CCD.crysol_stop_recursion = false;
    // reset sorting switches
    CCD.crysol_raw_sort_switch = 'RPSP';
    CCD.crysol_tagged_sort_switch = 'RPSP';
    CCD.crysol_cleaved_sort_switch = 'RPSP';

    hide_byclass('section3'); // it will be unhidden if polypeptides are selected
    hide_byclass('pp-crysol');
    hide_byclass('tag-crysol');
    hide_byclass('cleave-crysol');
    $('#crysol-save').css('display', 'none');

    $('.poly-check1').each(function (i) {
        if (this.checked) {
            selected_polypeptides_names.push(CCD.polypeptides[i][0]);
            selected_polypeptides_sequences.push(CCD.polypeptides[i][1]);
        }
    });
    $('.poly-check2').each(function (i) {
        if (this.checked) {
            tagged_selected_polypeptides_names.push(CCD.tagged_cleaved[i][0]);
            tagged_selected_polypeptides_sequences.push(CCD.tagged_cleaved[i][4]);
        }
    });
    $('.poly-check3').each(function (i) { // Difference from other polypeptides: can be cleaved by more than one proteases
        if (i > CCD.tagged_cleaved.length - 1) { // Tricky part: I cannot reset the i so I just subtract what would throw me out of bounds
            i = i - parseInt(i / CCD.tagged_cleaved.length) * CCD.tagged_cleaved.length; // could use mod function instead
        }
        if (this.checked) {
            cleaved_selected_polypeptides_names.push(CCD.tagged_cleaved[i][0] + ' with ' + this.name);
            $('.' + this.name + '-cleaved').each(function (j) {
                if (j === i) {
                    cleaved_selected_polypeptides_sequences.push($(this).html());
                }
            });
        }
    });
    if (selected_polypeptides_names.length === 0 &&
        tagged_selected_polypeptides_names.length === 0 &&
        cleaved_selected_polypeptides_names.length === 0) {
        $('#polypeptides-alerts').html('<div class="alert alert-danger"><a class="close" data-dismiss="alert">×</a>No polypeptides selected.</div>');
        return false;
    }

    // if pp selected, then show n request
    show_byclass('section3');

    if (selected_polypeptides_names.length !== 0) {
        CCD.crysol_raw = [];
        for (i = 0; i < selected_polypeptides_names.length; i += 1) {
            CCD.crysol_raw[i] = [selected_polypeptides_names[i], '', '', ''];
        }

        for (i = 0; i < 3; i += 1) { // pending = true for the raw predictions
            CCD.crysol_pending[0][i] = true;
        }

        rpsp_task(selected_polypeptides_sequences, CCD.crysol_raw, ' (raw)');
        proso2_task(selected_polypeptides_sequences, CCD.crysol_raw, ' (raw)');
        secret_task(selected_polypeptides_sequences, CCD.crysol_raw, ' (raw)');
    }
    if (tagged_selected_polypeptides_names.length !== 0 && $('#poly-tagged-toggle').attr('aria-pressed') === 'true') {
        CCD.crysol_tagged = [];
        for (i = 0; i < tagged_selected_polypeptides_names.length; i += 1) {
            CCD.crysol_tagged[i] = [tagged_selected_polypeptides_names[i], '', '', ''];
        }

        for (i = 0; i < 3; i += 1) { // pending = true for the tagged predictions
            CCD.crysol_pending[1][i] = true;
        }

        rpsp_task(tagged_selected_polypeptides_sequences, CCD.crysol_tagged, ' (tagged)');
        proso2_task(tagged_selected_polypeptides_sequences, CCD.crysol_tagged, ' (tagged)');
        secret_task(tagged_selected_polypeptides_sequences, CCD.crysol_tagged, ' (tagged)');
    }
    if (cleaved_selected_polypeptides_names.length !== 0 && $('#poly-cleaved-toggle').attr('aria-pressed') === 'true') {
        CCD.crysol_cleaved = [];
        for (i = 0; i < cleaved_selected_polypeptides_names.length; i += 1) {
            CCD.crysol_cleaved[i] = [cleaved_selected_polypeptides_names[i], '', '', ''];
        }

        for (i = 0; i < 3; i += 1) { // pending = true for the cleaved predictions
            CCD.crysol_pending[2][i] = true;
        }

        rpsp_task(cleaved_selected_polypeptides_sequences, CCD.crysol_cleaved, ' (cleaved)');
        proso2_task(cleaved_selected_polypeptides_sequences, CCD.crysol_cleaved, ' (cleaved)');
        secret_task(cleaved_selected_polypeptides_sequences, CCD.crysol_cleaved, ' (cleaved)');
    }

    $('#crysol-info').html('<div class="alert alert-info"><a class="close" data-dismiss="alert">×</a>Processing requests&nbsp;&nbsp;&nbsp;<img src="' + CCD.siteurl + 'static/img/bgLoad.gif' + '">');
}

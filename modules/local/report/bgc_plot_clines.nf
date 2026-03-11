process BGC_PLOT_CLINES {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(bgc_results)

    output:
        tuple val(meta), path("${meta.id}_gencline_plot_mqc.html"),         emit: cline_plot_html
        tuple val(meta), path("${meta.id}_gencline_param_scatter_mqc.html"), emit: scatter_plot_html
        tuple val(meta), path("${meta.id}_gencline_single_locus_mqc.html"),  emit: single_plot_html
        tuple val(meta), path("${meta.id}_cline_parameter_summary.tsv"), emit: summary_tsv

    script:
    def args = task.ext.args ?: ''

    """
    shopt -s nullglob

    center_matches=( "${bgc_results}"/*__gc_out__center.tsv )
    gradient_matches=( "${bgc_results}"/*__gc_out__gradient.tsv )
    locus_matches=( "${bgc_results}"/*__input_order__locus_order.txt )

    if (( \${#center_matches[@]} == 0 )); then
        echo "ERROR: Could not find cline center results file in ${bgc_results}" >&2
        echo "Contents of ${bgc_results}:" >&2
        ls -lah "${bgc_results}" >&2 || true
        exit 1
    fi

    if (( \${#gradient_matches[@]} == 0 )); then
        echo "ERROR: Could not find cline gradient results file in ${bgc_results}" >&2
        echo "Contents of ${bgc_results}:" >&2
        ls -lah "${bgc_results}" >&2 || true
        exit 1
    fi

    if (( \${#locus_matches[@]} == 0 )); then
        echo "ERROR: Could not find locus order file in ${bgc_results}" >&2
        echo "Contents of ${bgc_results}:" >&2
        ls -lah "${bgc_results}" >&2 || true
        exit 1
    fi

    CENTER_FILE="\${center_matches[0]}"
    GRADIENT_FILE="\${gradient_matches[0]}"
    LOCUS_FILE="\${locus_matches[0]}"

    echo "Using cline center file: \${CENTER_FILE}"
    echo "Using cline gradient file: \${GRADIENT_FILE}"
    echo "Using locus order file: \${LOCUS_FILE}"

    bgc_plot_clines.py \\
        --center            "\${CENTER_FILE}" \\
        --gradient          "\${GRADIENT_FILE}" \\
        --loci              "\${LOCUS_FILE}" \\
        --template-overlay  ${baseDir}/assets/multiqc_bgc_gencline_plot.html \\
        --template-scatter  ${baseDir}/assets/multiqc_bgc_gencline_scatter.html \\
        --template-single   ${baseDir}/assets/multiqc_bgc_gencline_single.html \\
        --out-overlay       "${meta.id}_gencline_plot_mqc.html" \\
        --out-scatter       "${meta.id}_gencline_param_scatter_mqc.html" \\
        --out-single        "${meta.id}_gencline_single_locus_mqc.html" \\
        --out-table ${meta.id}_cline_parameter_summary.tsv \\
        ${args}
    """
}

include { MULTIQC                } from '../../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_hybridclassification_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

include { PLOT_ADMIXTURE } from '../../modules/local/report/plot_admixture.nf'
include { NH_PLOT_CLASSIFICATIONS } from '../../modules/local/report/plot_newhybrids.nf'
include { NH_PLOT_TRACE } from '../../modules/local/report/newhybrids_trace.nf'
include { NH_SUMMARY_TABLE } from '../../modules/local/report/newhybrids_summary.nf'
include { NH_PLOT_SIMULATION } from '../../modules/local/report/newhybrids_simulation.nf'
include { NH_PLOT_SPATIAL } from '../../modules/local/report/plot_newhybrids_spatial.nf'
include { NH_PLOT_OUTLIERS } from '../../modules/local/report/plot_newhybrids_outliers.nf'
include { PLOT_TRIANGLE } from '../../modules/local/report/plot_triangle.nf'
include { BGC_MCMC_SUMMARY } from '../../modules/local/report/bgc_mcmc_summary.nf'
include { BGC_PLOT_HINDEX } from '../../modules/local/report/bgc_plot_hindex.nf'
include { BGC_PLOT_CLINES } from '../../modules/local/report/bgc_plot_clines.nf'

include { CUSTOMIZE_REPORT } from '../../modules/local/report/customize_report.nf'
include { STAGE_GEODATA_LAYERS } from '../../modules/local/report/stage_geodata.nf'

workflow GENERATE_REPORT {
    take:
    inds
    pops
    k2_clumpp
    bestK_clumpp
    nh_results
    nh_trace
    nh_map
    sim_results
    sim_trace
    sim_map
    triangle_popmap
    triangle_hindex
    triangle_hindex_fixed
    masked_samples
    candidate_map
    popmap
    speciesmap
    site_coords
    geo_data
    geo_data_dir
    bgc_output
    ch_versions

    main:

    ch_multiqc_files = Channel.empty()


    //Plot NH outliers based on Triangle test
    plot_mask_in = triangle_hindex
                    .join(triangle_popmap)
                    .join(nh_results)
                    .join(nh_map)
    NH_PLOT_OUTLIERS(
        plot_mask_in.map{ m, h, p, z, i -> tuple(m, h) },
        plot_mask_in.map{ m, h, p, z, i -> tuple(m, p) },
        plot_mask_in.map{ m, h, p, z, i -> tuple(m, z) },
        plot_mask_in.map{ m, h, p, z, i -> tuple(m, i) }
    )
    ch_multiqc_files = NH_PLOT_OUTLIERS.out.plot_html

    //Admixture barplots
    admix_inputs = k2_clumpp
                    .join(inds)
                    .join(pops)
                    .join(candidate_map)
    PLOT_ADMIXTURE(
        admix_inputs.map { m, k, i, p, c -> tuple(m, k) },
        admix_inputs.map { m, k, i, p, c -> tuple(m, i) },
        admix_inputs.map { m, k, i, p, c -> tuple(m, p) },
        admix_inputs.map { m, k, i, p, c -> tuple(m, c) }
    )
    ch_multiqc_files = ch_multiqc_files.join( PLOT_ADMIXTURE.out.admixture_html )
    ch_versions = ch_versions.mix( PLOT_ADMIXTURE.out.versions )

    //Triangle plot
    tri_input = triangle_hindex
                    .join( triangle_hindex_fixed )
                    .join( triangle_popmap )
                    .combine( popmap )
                    .map{ m, th, thf, tp, mp, p -> tuple(m, th, thf, tp, p) }
    PLOT_TRIANGLE(
        tri_input.map{ m, th, thf, tp, p -> tuple(m, th) },
        tri_input.map{ m, th, thf, tp, p -> tuple(m, thf) },
        tri_input.map{ m, th, thf, tp, p -> tuple(m, tp) },
        tri_input.map{ m, th, thf, tp, p -> tuple(m, p) }
    )
    ch_multiqc_files = ch_multiqc_files.join(PLOT_TRIANGLE.out.plot_html)

    // Power analysis plots
    sim_input = sim_results
                .join(sim_map)
    NH_PLOT_SIMULATION(
        sim_input.map{ m, si, sm -> tuple(m, si) },
        sim_input.map{ m, si, sm -> tuple(m, sm) }
    )
    ch_multiqc_files = ch_multiqc_files.join(NH_PLOT_SIMULATION.out.plot_html)
    ch_versions = ch_versions.mix( NH_PLOT_SIMULATION.out.versions )

    //NewHybrids pi trace plot
    NH_PLOT_TRACE(
        nh_trace
    )
    ch_multiqc_files = ch_multiqc_files.join(NH_PLOT_TRACE.out.plot_html)
    ch_versions = ch_versions.mix( NH_PLOT_TRACE.out.versions )

    //NewHybrids classifications
    nh_inputs = nh_results
                .join(nh_map)
                .join(masked_samples)
                .combine(popmap)
                .combine(speciesmap)
                .map{ m, nr, nm, ms, pm, p, sm, s -> tuple(m, nr, nm, ms, p, s) }
    NH_PLOT_CLASSIFICATIONS(
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, nr) },
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, nm) },
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, p) },
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, s) },
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, ms) },
    )
    ch_multiqc_files = ch_multiqc_files.join(NH_PLOT_CLASSIFICATIONS.out.plot_html)
    ch_versions = ch_versions.mix( NH_PLOT_CLASSIFICATIONS.out.versions )

    NH_SUMMARY_TABLE(
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, nr) },
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, nm) },
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, p) },
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, s) },
        nh_inputs.map{m, nr, nm, ms, p, s -> tuple(m, ms) }
    )
    ch_multiqc_files = ch_multiqc_files.join(NH_SUMMARY_TABLE.out.table_json)
    ch_multiqc_files = ch_multiqc_files.join(NH_SUMMARY_TABLE.out.table_json_masked)
    ch_versions = ch_versions.mix( NH_SUMMARY_TABLE.out.versions )

    //
    // Map of newhybrids results
    //
    if (params.site_coords){
        if (params.geo_data_config){

            STAGE_GEODATA_LAYERS( geo_data, geo_data_dir )

            ch_spatial_inputs = nh_results
                                .join( nh_map )
                                .join( masked_samples )
                                .combine( popmap )
                                .combine( site_coords )
                                .combine( STAGE_GEODATA_LAYERS.out.geo_data_dir )
                                .map{
                                    m, nr, nm, mask, mp, p, ms, s, mg, g -> tuple( m, nr, nm, mask, p, s, g )
                                }
            NH_PLOT_SPATIAL(
                ch_spatial_inputs.map{m, nr, nm, mask, p, s, g -> tuple(m, nr) },
                ch_spatial_inputs.map{m, nr, nm, mask, p, s, g -> tuple(m, nm) },
                ch_spatial_inputs.map{m, nr, nm, mask, p, s, g -> tuple(m, p) },
                ch_spatial_inputs.map{m, nr, nm, mask, p, s, g -> tuple(m, s) },
                ch_spatial_inputs.map{m, nr, nm, mask, p, s, g -> tuple(m, g) },
                ch_spatial_inputs.map{m, nr, nm, mask, p, s, g -> tuple(m, mask) }
            )
            ch_multiqc_files = ch_multiqc_files.join(NH_PLOT_SPATIAL.out.plot_html)
            ch_versions = ch_versions.mix( NH_PLOT_SPATIAL.out.versions )
        } else {
            ch_spatial_inputs = nh_results
                .join( nh_map )
                .join( masked_samples )
                .combine( popmap )
                .combine( site_coords )
                .map { m, nr, nm, mask, mp, p, ms, s -> tuple( m, nr, nm, mask, p, s) }
            NH_PLOT_SPATIAL(
                ch_spatial_inputs.map { m, nr, nm, mask,  p, s -> tuple(m, nr) },
                ch_spatial_inputs.map { m, nr, nm, mask,  p, s -> tuple(m, nm) },
                ch_spatial_inputs.map { m, nr, nm, mask,  p, s -> tuple(m, p) },
                ch_spatial_inputs.map { m, nr, nm, mask,  p, s -> tuple(m, s) },
                tuple( [], [] ),
                ch_spatial_inputs.map { m, nr, nm, mask,  p, s -> tuple(m, mask) }
            )
            ch_multiqc_files = ch_multiqc_files.join(NH_PLOT_SPATIAL.out.plot_html)
            ch_versions = ch_versions.mix( NH_PLOT_SPATIAL.out.versions )
        }
    }

    //
    // Optional genomic clines handling
    //
    if (params.run_bgc){

        // MCMC trace plot
        // Not planned for now (inflates report size too much)

        // MCMC summary stats
        BGC_MCMC_SUMMARY( bgc_output )
        ch_multiqc_files = ch_multiqc_files.join(BGC_MCMC_SUMMARY.out.table_json)

        // Hybrid indices
        BGC_PLOT_HINDEX( bgc_output )
        ch_multiqc_files = ch_multiqc_files.join(BGC_PLOT_HINDEX.out.plot_html)

        // Genomic clines and alpha-beta plot
        BGC_PLOT_CLINES( bgc_output )
        ch_multiqc_files = ch_multiqc_files.join(BGC_PLOT_CLINES.out.cline_plot_html)
        ch_multiqc_files = ch_multiqc_files.join(BGC_PLOT_CLINES.out.scatter_plot_html)
        ch_multiqc_files = ch_multiqc_files.join(BGC_PLOT_CLINES.out.single_plot_html)
    }


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // Prepare MultiQC files
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    //
    // Merge MultiQC files with report plots
    //
    ch_multiqc_files = ch_multiqc_files
        .combine(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .combine(ch_collated_versions)
        .combine(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
    )

    MULTIQC (
        ch_multiqc_files,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    CUSTOMIZE_REPORT(
        MULTIQC.out.report
    )

    emit:
    multiqc_report  = MULTIQC.out.report.toList()
    versions        = ch_collated_versions
}

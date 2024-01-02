"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""
with open(config["antibody_escape_config"]) as f:
    antibody_escape_config = yaml.safe_load(f)

antibody_escape_pdbs = antibody_escape_config["antibody_escape_PDBs"]

Env_chains_by_pdb = antibody_escape_config["env_chains_by_PDB"]

chains_to_exclude = antibody_escape_config["chains_to_exclude"]

antibody_gv = [x for x in antibody_escape_config['avg_antibody_escape']]

rule spatial_distances:
    """Get spatial distances from PDB."""
    input: 
        pdb="data/env_trimer_6UDJ.pdb",
    output:
        csv="results/spatial_distances/6udj.csv",
    params:
        target_chains=["A", "B", "C"],
    log:
        log="results/logs/spatial_distances.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    script:
        "scripts/spatial_distances.py"

rule color_PDB_structures:
    """Assign b factor values to PDB structures based on escape"""
    input: 
        average_escape_model = "results/antibody_escape/averages/{antibody}_mut_effect.csv",
        input_pdb_file = lambda wc: os.path.join(
            antibody_escape_config["PDB_structures"], 
            antibody_escape_pdbs[wc.antibody]
        ),
    output:
        output_pdb_file_name = os.path.join("results/escape_PDBs/", "{antibody}" + ".pdb"),
    params: 
        env_chains = lambda wc: Env_chains_by_pdb[antibody_escape_pdbs[wc.antibody]],
        output_file_sub_name = os.path.join("results/escape_PDBs/", "{antibody}"),
    log:
        os.path.join("results/logs/", "antibody_escape_pdbs_{antibody}.txt"),
    script:
        "scripts/color_pdb_structures.py"

rule format_dms_viz:
    """Format the data for input into dms-viz"""
    input:
        average_escape_model = rules.avg_escape.output.effect_csv,
        input_pdb_file = lambda wc: os.path.join(
            antibody_escape_config["PDB_structures"], 
            antibody_escape_pdbs[wc.antibody]
        ),
        site_map="data/site_numbering_map.csv",
    output: 
        output_json_file_name = os.path.join("results/dms-viz/", "{assay}", "{antibody}" + ".json"),
    params:
        env_chains = lambda wc: Env_chains_by_pdb[antibody_escape_pdbs[wc.antibody]],
        exclude_chains = lambda wc: chains_to_exclude[wc.antibody],
        name="{antibody}",
    log:
        os.path.join("results/logs/", "dms-viz_file_{assay}_{antibody}.txt"),
    conda:
        "dms-viz.yml"
    shell:
        """
        configure-dms-viz format \
            --name {params.name} \
            --input {input.average_escape_model} \
            --metric  "escape_mean" \
            --metric-name "Escape" \
            --condition "epitope" \
            --condition-name "Epitope" \
            --structure {input.input_pdb_file} \
            --sitemap {input.site_map} \
            --output {output.output_json_file_name} \
            --included-chains "{params.env_chains}" \
            --excluded-chains "{params.exclude_chains}" \
            --tooltip-cols "{{'times_seen': '# Obsv'}}" \
            --filter-cols "{{'times_seen': 'Times Seen'}}" \
            &> {log}
        """

rule correlate_escape:
    """Correlate antibody escape"""
    input: 
        "results/antibody_escape/averages/1-18_mut_effect.csv",
        "results/antibody_escape/averages/3BNC117_mut_effect.csv",
        "results/antibody_escape/averages/04-A06_mut_effect.csv",
        "results/antibody_escape/averages/N6_LS_mut_effect.csv",
        "results/antibody_escape/averages/VRC07_523_LS_mut_effect.csv",
        nb="notebooks/escape_correlations.ipynb",
    output:
        nb="results/notebooks/escape_correlations.ipynb"
    log:
        "results/logs/escape_correlations.txt"
    shell:
        """
        papermill {input.nb} {output.nb} \
            &> {log}
        """

# Files (Jupyter notebooks, HTML plots, or CSVs) that you want included in
# the HTML docs should be added to the nested dict `docs`:
docs["Additional analysis-specific files"] = {
    "Reference to sequential site-numbering map": config["site_numbering_map"],
    "Correlations of antibody escape effects": "results/notebooks/escape_correlations.ipynb",
}

# If you want to make other output files from your rules target files for
# the pipeline, add them to `other_target_files` list
for file in expand(rules.format_dms_viz.output.output_json_file_name, antibody=antibody_gv, assay='antibody_escape'):
    other_target_files.append(file)
for file in expand(rules.color_PDB_structures.output.output_pdb_file_name, antibody=antibody_gv, assay='antibody_escape'):
    other_target_files.append(file)


def get_final_output():
    final_output = expand(
        "results/diffexp/{contrast}.diffexp.symbol.tsv",
        contrast=config["diffexp"]["contrasts"],
    )
    final_output.append("results/deseq2/normcounts.symbol.tsv")
    final_output.append("results/counts/all.symbol.tsv")
    final_output.append("results/qc/multiqc_report.html")

    if config["pca"]["activate"]:
        # get all the variables to plot a PCA for
        pca_variables = list(config["diffexp"]["variables_of_interest"])
        if config["diffexp"]["batch_effects"]:
            pca_variables.extend(config["diffexp"]["batch_effects"])
        if config["pca"]["labels"]:
            pca_variables.extend(config["pca"]["labels"])
        final_output.extend(
            expand("results/pca.{variable}.svg", variable=pca_variables)
        )
    return final_output

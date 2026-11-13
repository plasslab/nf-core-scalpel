# plotExpression

Plot Relative Isoform or Feature Expression Across Groups

```{raw} html



<h3>Description</h3>

<p>Generates boxplots of relative expression for selected features (e.g., isoforms or genes) across groups in a Seurat object.
</p>


<h3>Usage</h3>

<pre><code class='language-R'>plotExpression(
  seurat_obj,
  features,
  group_var,
  plot = "group",
  assay = "RNA",
  layer = "counts",
  log_normalization = TRUE,
  only_expressed_cells = FALSE
)
</code></pre>


<h3>Arguments</h3>

<table>
<tr><td><code id="seurat_obj">seurat_obj</code></td>
<td>
<p>A Seurat object containing single-cell RNA-seq data.</p>
</td></tr>
<tr><td><code id="features">features</code></td>
<td>
<p>Character vector of features (e.g., genes or isoforms) to plot.</p>
</td></tr>
<tr><td><code id="group_var">group_var</code></td>
<td>
<p>The metadata column in the Seurat object to group cells by.</p>
</td></tr>
<tr><td><code id="plot">plot</code></td>
<td>
<p>Character. Whether to plot by group (<code>"group"</code>) or by features (<code>"features"</code>). Default is <code>"group"</code>.</p>
</td></tr>
<tr><td><code id="assay">assay</code></td>
<td>
<p>Character. The Seurat assay to use. Default is <code>"RNA"</code>.</p>
</td></tr>
<tr><td><code id="layer">layer</code></td>
<td>
<p>Character. The layer to use (layers split per sample are joined). Default is <code>"counts"</code>.</p>
</td></tr>
<tr><td><code id="log_normalization">log_normalization</code></td>
<td>
<p>Logical. Whether to log-normalize the expression values (<code>TRUE</code>) or plot raw counts (<code>FALSE</code>). Default is <code>TRUE</code>.</p>
</td></tr>
<tr><td><code id="only_expressed_cells">only_expressed_cells</code></td>
<td>
<p>Logical. Plot only the cells expressing the feature. Default is <code>FALSE</code>.</p>
</td></tr>
</table>


<h3>Details</h3>

<p>The function extracts counts for the specified features from the specified assay in the Seurat object. It plots boxplots of expression values, grouped either by metadata group or by features, and can display either raw or log-normalized counts.
</p>


<h3>Value</h3>

<p>A ggplot2 object showing boxplots of expression values for the selected features across the specified groups.
</p>


<h3>Examples</h3>

<pre><code class='language-R'>## Not run: 
plotExpression(
  seurat_obj = my_seurat_obj,
  features = c("Isoform1", "Isoform2"),
  group_var = "Cluster",
  plot = "group",
  assay = "RNA",
  log_normalization = TRUE
)

## End(Not run)

</code></pre>

</main>

```

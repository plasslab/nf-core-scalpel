# weightedUTR_cell

Calculate Weighted UTR Length per Cell

```{raw} html



<h3>Description</h3>

<p>This function computes the weighted 3' UTR length per cell using the isoform quantification of a Seurat object
and a gene annotation with UTR and CDS information (3' UTR length of each isoform).
</p>


<h3>Usage</h3>

<pre><code class='language-R'>weightedUTR_cell(
  seurat_obj,
  annotation_gr,
  group.by,
  features = NULL,
  assay = "RNA",
  method = "centered",
  min.counts = NULL
)
</code></pre>


<h3>Arguments</h3>

<table>
<tr><td><code id="seurat_obj">seurat_obj</code></td>
<td>
<p>A Seurat object containing single-cell RNA-seq data.</p>
</td></tr>
<tr><td><code id="annotation_gr">annotation_gr</code></td>
<td>
<p>A GRanges object containing gene annotation, including UTR and CDS regions.</p>
</td></tr>
<tr><td><code id="group.by">group.by</code></td>
<td>
<p>The metadata column in the Seurat object to group cells by.</p>
</td></tr>
<tr><td><code id="features">features</code></td>
<td>
<p>A vector of features (genes/transcripts) to subset the Seurat object by. Default is NULL.</p>
</td></tr>
<tr><td><code id="assay">assay</code></td>
<td>
<p>The assay to use for expression values, default is &quot;RNA&quot;.</p>
</td></tr>
<tr><td><code id="method">method</code></td>
<td>
<p><code>"centered"</code> (default) or <code>"length"</code>, see Details.</p>
</td></tr>
<tr><td><code id="min.counts">min.counts</code></td>
<td>
<p>Isoforms with <code>min.counts</code> counts or less in a cell are not used in this cell. Default:
0 with <code>method = "centered"</code>, 3 with <code>method = "length"</code>.</p>
</td></tr>
</table>


<h3>Details</h3>

<p>For each gene and cell, the weighted 3' UTR length is the sum over its isoforms of the isoform fraction of the gene
counts x the isoform 3' UTR length.
</p>
<p><code>method = "centered"</code>: only genes whose isoforms differ in 3' UTR length are used, and the weighted 3' UTR length
of each gene in the cell is centered on the weighted 3' UTR length of the gene in all the cells of the object;
<code>wUTRlength.cell</code> is the mean of these differences over the genes of the cell (nt; &gt; 0: longer 3' UTRs than the
average cell, &lt; 0: shorter). It does not depend on which genes are detected in the cell, i.e. on the sequencing depth.
</p>
<p><code>method = "length"</code>: <code>wUTRlength.cell</code> is the mean over the genes of the cell of their weighted 3' UTR
length (nt). Genes with a single 3' UTR length contribute their fixed length: the value depends on the genes detected
in the cell, which increase with the sequencing depth.
</p>
<p>Counts split in several layers (one per sample) are joined; computed on sparse matrices.
</p>


<h3>Value</h3>

<p>A data.table with one row per cell: <code>cells</code>, <code>group</code>, <code>wUTRlength.cell</code> and
<code>n_genes</code> (genes used in the cell).
</p>


<h3>Examples</h3>

<pre><code class='language-R'>## Not run: 
weightedUTR_cell(seurat_obj = my_seurat_obj, annotation_gr = my_annotation, group.by = "Cluster")

## End(Not run)

</code></pre>

</main>

```

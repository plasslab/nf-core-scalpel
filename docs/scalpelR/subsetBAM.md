# subsetBAM

Subset BAM File by Barcodes

```{raw} html



<h3>Description</h3>

<p>Subsets a BAM file to include only reads with specified barcode tags, writing the result to a new BAM file.
</p>


<h3>Usage</h3>

<pre><code class='language-R'>subsetBAM(
  bam_file,
  output_file,
  barcodes,
  tag,
  samtools_bin = "samtools",
  cores = 4
)
</code></pre>


<h3>Arguments</h3>

<table>
<tr><td><code id="bam_file">bam_file</code></td>
<td>
<p>Path to the input BAM file to be subset.</p>
</td></tr>
<tr><td><code id="output_file">output_file</code></td>
<td>
<p>Path to the output BAM file to be written.</p>
</td></tr>
<tr><td><code id="barcodes">barcodes</code></td>
<td>
<p>Character vector of cell barcodes to retain, as in the BAM tag (SCALPEL iDGE Seurat object:
<code>Barcodes</code> column; cell names are <code>&lt;sample&gt;_&lt;barcode&gt;</code>).</p>
</td></tr>
<tr><td><code id="tag">tag</code></td>
<td>
<p>Character. The tag in the BAM file identifying the barcode (e.g., <code>"BC"</code>).</p>
</td></tr>
<tr><td><code id="samtools_bin">samtools_bin</code></td>
<td>
<p>Character. Path to the <code>samtools</code> binary. Default is <code>"samtools"</code>.</p>
</td></tr>
<tr><td><code id="cores">cores</code></td>
<td>
<p>Integer. Number of cores to use for filtering. Default is <code>4</code>.</p>
</td></tr>
</table>


<h3>Details</h3>

<p>The function writes the selected barcodes to a temporary file and uses <code>samtools view</code> to extract reads with matching barcode tags into a new BAM file, followed by indexing.
</p>


<h3>Value</h3>

<p>NULL. The function is called for its side effects (writing a subset BAM and its index file).
</p>


<h3>Examples</h3>

<pre><code class='language-R'>## Not run: 
subsetBAM(
  bam_file = "input.bam",
  output_file = "subset.bam",
  barcodes = c("AAACCTGAGC", "AAACCTGCTT"),
  tag = "BC",
  samtools_bin = "samtools",
  cores = 4
)

## End(Not run)

</code></pre>

</main>

```

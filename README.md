<h1 style="color:#4da6ff;">Vina-SBVS-LBVS-ADME: Full Workflow Documentation</h1>

<p>A complete protocol for protein preparation, ligand preparation, pocket detection, AutoDock Vina docking, post-processing, ADME integration, ligand efficiency scoring, and visualization.</p>

<hr>

<h2 style="color:#4da6ff;">📌 Table of Contents</h2>

<ul>
  <li><a href="#overview">Overview</a></li>
  <li><a href="#system-requirements">System Requirements</a></li>
  <li><a href="#software-installation">Software Installation</a></li>
  <li><a href="#recommended-directory-structure">Recommended Directory Structure</a></li>
  <li><a href="#protein-preparation">Protein Preparation</a></li>
  <li><a href="#ligand-library-preparation">Ligand Library Preparation</a></li>
  <li><a href="#ligand-protonation-3d-generation--minimization">Ligand Protonation, 3D Generation & Minimization</a></li>
  <li><a href="#adme-descriptor-acquisition">ADME Descriptor Acquisition</a></li>
  <li><a href="#running-the-main-bash-workflow">Running the Main Bash Workflow</a></li>
  <li><a href="#pocket-detection--docking-workflow">Pocket Detection & Docking Workflow</a></li>
  <li><a href="#r-based-post-processing-workflow">R-Based Post-Processing Workflow</a></li>
  <li><a href="#thermodynamics--ligand-efficiency-metrics">Thermodynamics & Ligand Efficiency Metrics</a></li>
  <li><a href="#adme-drug-likeness-rules--boiled-egg">ADME, Drug-Likeness Rules & BOILED-Egg</a></li>
  <li><a href="#visualization-outputs">Visualization Outputs</a></li>
  <li><a href="#troubleshooting">Troubleshooting</a></li>
  <li><a href="#references">References</a></li>
</ul>

<hr>

<h2 id="overview" style="color:#4da6ff;">Overview</h2>

<p>This repository provides a unified, fully automated computational workflow that includes:</p>

<ul>
  <li>Viral protein preparation</li>
  <li>Ligand preparation and format conversion</li>
  <li>Pocket detection using Concavity</li>
  <li>Structure-based virtual screening via AutoDock Vina</li>
  <li>Automated log parsing</li>
  <li>Binding free energy (ΔG) extraction</li>
  <li>ΔG → Kd thermodynamic conversion</li>
  <li>ADME descriptor integration</li>
  <li>Ligand efficiency metrics (LE, LLE, FQ)</li>
  <li>Heatmaps, 3D affinity surfaces, and PCA-ready datasets</li>
</ul>

<hr>

<h2 id="system-requirements" style="color:#4da6ff;">System Requirements</h2>

<ul>
  <li>Linux (Ubuntu recommended)</li>
  <li>Bash shell</li>
  <li>AutoDock Vina</li>
  <li>MGLTools</li>
  <li>OpenBabel</li>
  <li>AmberTools / Antechamber</li>
  <li>Concavity tool</li>
  <li>Python 3</li>
  <li>R + RStudio</li>
</ul>

<hr>

<h2 id="software-installation" style="color:#4da6ff;">Software Installation</h2>

<p>Install the following tools:</p>

<ul>
  <li>AutoDock Vina – https://github.com/ccsb-scripps/AutoDock-Vina/releases</li>
  <li>MGLTools – https://autodock.scripps.edu</li>
  <li>OpenBabel – https://openbabel.org</li>
  <li>AmberTools – https://ambermd.org</li>
  <li>R – https://cran.r-project.org</li>
  <li>RStudio – https://posit.co</li>
  <li>APBS (optional) – https://www.poissonboltzmann.org</li>
</ul>

<p><b>Note:</b> Ensure MGLTools <i>Utilities24</i> path is correctly set in <code>VINA_with_conversion.sh</code>.</p>

<hr>

<h2 id="recommended-directory-structure" style="color:#4da6ff;">Recommended Directory Structure</h2>

<pre><code>~/Vina_Workspace/
    receptors/
        WNV_E/
            receptor_WNV_E.pdb
        YFV_NS5/
            receptor_YFV_NS5.pdb
        ...
    compounds/
        ligand1.sdf
        ligand2.sdf
        ...
</code></pre>

<p>Final outputs (heatmaps, CSVs, Kd tables, etc.) are generated automatically.</p>

<hr>

<h2 id="protein-preparation" style="color:#4da6ff;">Protein Preparation</h2>

<p>Create directories:</p>

<pre><code>mkdir ~/Vina_Workspace/receptors/WNV_E
</code></pre>

<p>Place receptor:</p>

<pre><code>receptor_WNV_E.pdb
</code></pre>

<p>Pipeline automatically:</p>
<ul>
  <li>Identifies target</li>
  <li>Generates PDBQT</li>
  <li>Detects pockets</li>
</ul>

<hr>

<h2 id="ligand-library-preparation" style="color:#4da6ff;">Ligand Library Preparation</h2># Automated-Virtual-Screening-and-ADME-Workflow-Vina-R

<p>Place all ligands in:</p>

<pre><code>~/Vina_Workspace/compounds/
</code></pre>

<p>Ligands lacking 3D coordinates should be renamed with prefix <b>3D_</b>.</p>

<hr>

<h2 id="ligand-protonation-3d-generation--minimization" style="color:#4da6ff;">Ligand Protonation, 3D Generation & Minimization</h2>

<h3 style="color:#33bbff;">pH adjustment</h3>

<pre><code>cd ~/Vina_Workspace/compounds
for f in *.sdf; do
    tmp="${f}.tmp"
    obabel -i sdf "$f" -o sdf -O "$tmp" -p 7
    mv "$tmp" "$f"
done
</code></pre>

<h3 style="color:#33bbff;">3D Generation + Minimization</h3>

<pre><code>cd ~/Vina_Workspace/compounds
for f in *.sdf; do
    tmp="${f}.tmp"
    obabel -i sdf "$f" -o sdf -O "$tmp" \
        --conformer --nconf 1000 --weighted \
        --minimize --steps 1000
    mv "$tmp" "$f"
done
</code></pre>

<h3 style="color:#33bbff;">SDF → MOL2 Conversion</h3>

<pre><code>cd ~/Vina_Workspace/compounds
for file in *.sdf; do
    antechamber -i "$file" -fi sdf \
        -o "${file%.sdf}_antechamber.mol2" -fo mol2 \
        -c bcc -at gaff2
done
</code></pre>

<p><b style="color:red;">IMPORTANT:</b> Ligand names must NOT contain underscores.</p>

<hr>

<h2 id="adme-descriptor-acquisition" style="color:#4da6ff;">ADME Descriptor Acquisition</h2>

<p>Convert MOL2 to SMILES:</p>

<pre><code>cd ~/Vina_Workspace/compounds
obabel *_antechamber.mol2 -O adme_input.smiles -osmi --title
</code></pre>

<p>Upload to SwissADME → export CSV:</p>

<pre><code>adme_features_for_efficiency.csv
</code></pre>

<hr>

<h2 id="running-the-main-bash-workflow" style="color:#4da6ff;">Running the Main Bash Workflow</h2>

<p>Enable script:</p>

<pre><code>chmod +x VINA_with_conversion.sh
</code></pre>

<p>Run:</p>

<pre><code>./VINA_with_conversion.sh
</code></pre>

<p>Pipeline performs:</p>

<ul>
  <li>Ligand & receptor conversion</li>
  <li>Pocket detection</li>
  <li>Grid box creation</li>
  <li>Docking</li>
  <li>Log organization</li>
  <li>Report generation</li>
</ul>

<hr>

<h2 id="pocket-detection--docking-workflow" style="color:#4da6ff;">Pocket Detection & Docking Workflow</h2>

<ul>
  <li>prepare_receptor4.py</li>
  <li>prepare_ligand4.py</li>
  <li>Concavity analysis</li>
  <li>Pocket centroid calculation</li>
  <li>Configuration generation</li>
</ul>

<p>Docking params:</p>

<pre><code>num_modes = 9
exhaustiveness = 8
energy_range = 4
</code></pre>

<hr>

<h2 id="r-based-post-processing-workflow" style="color:#4da6ff;">R-Based Post-Processing Workflow</h2>

<p>Install libs:</p>

<pre><code>install.packages(c(
"dplyr","tidyr","stringr","ComplexHeatmap",
"circlize","plotly","ggplot2","viridis",
"rcdk","janitor","htmlwidgets"
))
</code></pre>

<p>Run:</p>

<pre><code>VINA_WORKSPACE <- path.expand("~/Vina_Workspace")
source("Main_analysis.R")
</code></pre>

<p>Generates matrices, heatmaps, 3D plots, rankings and CSVs.</p>

<hr>

<h2 id="thermodynamics--ligand-efficiency-metrics" style="color:#4da6ff;">Thermodynamics & Ligand Efficiency Metrics</h2>

<p>ΔG → Kd:</p>

<pre><code>Kd = exp(ΔG / (R * T))
T = 310.15 K
R = 1.987e-3
</code></pre>

<p>LE, LLE, FQ formulas included.</p>

<hr>

<h2 id="adme-drug-likeness-rules--boiled-egg" style="color:#4da6ff;">ADME, Drug-Likeness Rules & BOILED-Egg</h2>

<p>Evaluates Lipinski, Ghose, Egan, Veber + BOILED-Egg (HIA/BBB).</p>

<hr>

<h2 id="visualization-outputs" style="color:#4da6ff;">Visualization Outputs</h2>

<ul>
  <li>ΔG heatmaps</li>
  <li>3D affinity landscapes</li>
  <li>LLE × pKd scatterplots</li>
  <li>BOILED-Egg diagrams</li>
</ul>

<hr>

<h2 id="troubleshooting" style="color:#4da6ff;">Troubleshooting</h2>

<ul>
  <li>Underscores in ligand names → remove</li>
  <li>MGLTools path wrong → fix Utilities24</li>
  <li>Pockets empty → check receptor integrity</li>
  <li>ADME merge fails → standardize names</li>
</ul>

<hr>

<h2 id="references" style="color:#4da6ff;">References</h2>

<ul>
  <li>AutoDock Vina</li>
  <li>SwissADME</li>
  <li>Ligand Efficiency metrics</li>
  <li>Drug-likeness literature</li>
</ul>

<hr>

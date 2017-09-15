# EecSeq Lab Protocol

The Expressed Exome Capture Sequencing protocol is designed to create exome capture probes directly from RNA.  The probes are then used from hybrid capture of exome DNA sequences, allowing for genotyping of alleles at expressed genes.

# Outline

[RNA Prep](#rna-prep)
* [RNA Extraction](#extract-rna-from-individuals-to-be-used-for-probes)
* [Quantify](#quantify-all-rna-samples)
* [Visualize](#visualize-rna-on-bioanalyzer)
* [mRNA Library Prep](#stranded-mrna-seq-library-prep)
	* [Anneal adapters](#anneal-rna-adapters)
	* [mRNA capture](#mrna-capture)
	* [mRNA Elution, Fragmentation, and Priming](#mrna-elution,-fragmentation-and-priming)
	* [1st Strand Synthesis](#1st-strand-synthesis)
	* [2nd Strang Synthesis](#2nd-strand-synthesis-and-marking)
	* [A-tailing](#a-tailing)
	* [Adapter Ligation](#adapter-ligation)
	* [Library Enrichment](#library-amplificiation)
* [Quant](#-quant-libraries)
* [DSN Normalization](#dsn-normalization)
	* [Enrichment](#pcr-enrichment)
* [Quant Libraries](#quant-libraries)
* [Split libraries](#split-finished-cdna-library-for-each-sample-into-two-vials-8-vials-total)

[Probe Synthesis](#probe-synthesis)
* [Remove sequencing adapters](#remove-adapters-from-cdna)
* [Remove 5' and 3' overhangs](#remove-5'-and-3'-overhangs)
* [Biotin Labeling](#biotin-labeling)

[Genomic Library Prep](#preparation-of-whole-genome-libraries-using-kapa-hyperplus-kit)
* [Anneal Adapters](#anneal-adapters)
* [End Repair](#end-repair)
* [Adapter ligation](#adapter-ligation)
* [Quantification](#quant-samples)
* [Amplification](#library-amplification)

[Hybridization and Capture](#hybridization-and-capture)
* [Hybridization](#hybridization)
	* [Preparation of DynaBeads](#preparation-of-dynabeads)
	* [Washes](#washes)
* [Library Re-amplification](#library-re-amplification)
* [Quant Final Libraries](#quant-samples)
* [Verify](#verify)



---



## RNA Prep

### Extract RNA from individuals to be used for probes
*Refer to manual during procedure (steps below are for notes and comments)*
### Using unmodified TRI Reagent Protocol [LINK](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/9738M_D.pdf)(Below are summary steps)

#### Reagents and supplies
* Nuclease-free Water
* 1-bromo-3-chloropropane (BCP; recommended, e.g., MRC, Cat #BP 151), or chloroform without added isoamyl alcohol 100% 
* isopropanol, ACS grade or better
* 100% ethanol, ACS grade or better

#### Equipment
* Tissue homogenizer with replacable probes (such as Qiagen TissueRuptor®)
* Appropriately sized RNase-free centrifuge tubes with secure closures, compatible with phenol/chloroform (polypropyl- ene, or polyallomer), and capable of withstanding centrifugal forces of 12,000G
* Centrifuge capable of 12,000 x g

**Notes before starting**
* Place all pipettes, tips, and supplies inside the hood and expose to UV for 30 minutes. (Do not put your samples in the UV hood!)
* Prepare 75% ethanol by mixing 250 μL nuclease-free water with 750 μL 100% ethanol per mL of TRI Reagent solution to be used. Include 10% overage to ensure a sufficient volume.

#### Procedure
* Homogenize tissue samples in 10–20 volumes TRI Reagent solution. Homogenize cultured cells in 1 mL TRI Reagent solution per 5–10 x 106 cells, or per 10 cm2 culture dish area.
* Incubate the homogenate for 5 min at room temp.
* (Optional) Centrifuge at 12,000 xg for 10 min at 4°C and transfer the supernatant to a fresh tube.
* Add 100 μL BCP per 1 mL of TRI Reagent solution, mix well, and incubate at room temp for 5–15 min.
	* Add 100 μL BCP per 1 mL of TRI Reagent solution used for homogenization. Alternatively, 200 μL of chloroform (without isoamyl alcohol) can be used in place of BCP.
	* Cap the tubes tightly and shake vigorously for 15 sec.
	* Incubate the mixture at room temperature for 5–15 min.
* Centrifuge at 12,000 x g for 10–15 min at 4°C, then transfer the aqueous phase to a fresh tube.
	* Centrifuge at 12,000 x g for 10–15 min at 4°C. Centrifugation at temperatures > 8°C may cause some DNA to partition in the aqueous phase.
	* Transfer the aqueous phase (colorless top layer) to a fresh tube.
* Add 500 μL of isopropanol per 1 mL of TRI Reagent solution, vortex for 5–10 sec, and incubate at room temp for 5–10 min.
	* Add 500 μL of isopropanol per 1 mL of TRI Reagent solution used for sample homogenization.
	* Vortex at moderate speed for 5–10 sec.
	* Incubate the samples at room temp for 5–10 min.
* Centrifuge at 12,000 x g for 8 min at 4–25°C, and discard the supernatant.
	* Centrifuge at 12,000 x g for 8 min at 4–25°C.
	* Carefully remove the supernatant without disturbing the pellet.
		* Precipitated RNA forms a gel-like or white pellet on the side and bottom of the tube.
* Add 1 mL of 75% ethanol per 1 mL of TRI Reagent solution
* Centrifuge at 7,500 x g for 5 min, remove the ethanol, and briefly air dry the RNA pellet.
	* Centrifuge at 7,500 x g for 5 min at 4–25°C.  
		* If the precipitated RNA floats or does not form a compact pellet, repeat the centrifugation at 12,000 x g for 5 min to consolidate the pellet at the bottom of the tube.
	* Remove the ethanol wash without disturbing the pellet.
	* Remove all residual ethanol by centrifuging again briefly and removing the ethanol that collects with a fine tip pipette. **Complete removal of ethanol is necessary for the RNA to perform well in downstream applications.**
	* Air dry the RNA pellet for 3–5 min.
		* Do not completely dry the RNA pellet as this will greatly decrease its solubility. 
* Dissolve RNA in the buffer of your choice.
	* Dissolve RNA in THE RNA Storage Solution (P/N AM7000, AM7001), Nuclease-free Water, or your choice of buffer‡ by passing the solution a few times through a pipette tip or by vigorous vortexing
		* The resuspension volume is determined by the size of the RNA pellet. 3–5 mm pellets typically require 300–500 μL. If necessary, increase the resuspension volume or incubate at 55–60°C to completely dissolve the pellet.
	* Store at 4°C for immediate analysis. For long-term storage, store at –70°C or colder.
	
### Quantify all RNA samples
Results will be used for calibration points during library generation
Refer to manual during procedure (steps below are for notes and comments)

#### Reagents and supplies
* Qubit® RNA HS Assay Kit (ThermoFisher Q32852)
* Microcentrifuge tubes for florescence (Fisher Catalog # 07-200-183)

#### Equipment
* Qubit® 3.0 Flourometer

#### Procedure (Standard HS RNA protocol)
* Set up the required number of 0.5-mL tubes for standards and samples. The Qubit® RNA HS Assay requires 2 standards.
* Label the tube lids.
* Prepare the Qubit® working solution by diluting the Qubit® RNA HS Reagent 1:200 in Qubit® RNA HS Buffer. Use a clean plastic tube each time you prepare Qubit® working solution. **Do not mix the working solution in a glass container.**
* Add 190 μL of Qubit® working solution to each of the tubes used for standards.
* Add 10 μL of each Qubit® standard to the appropriate tube, then mix by vortexing 2–3 seconds. Be careful not to create bubbles.
* Add Qubit® working solution to individual assay tubes so that the final volume in each tube after adding sample is 200 μL.
* Add each sample to the assay tubes containing the correct volume of Qubit® working solution, then mix by vortexing 2–3 seconds. The final volume in each tube should be 200 μL.
* Allow all tubes to incubate at room temperature for 2 minutes.
* On the Home screen of the Qubit® 3.0 Fluorometer, press RNA, then select RNA: High Sensitivity as the assay type. The “Read standards” screen is displayed. Press Read Standards to proceed.
* Insert the tube containing Standard #1 into the sample chamber, close the lid, then press Read standard. When the reading is complete (~3 seconds), remove Standard #1.
* Insert the tube containing Standard #2 into the sample chamber, close the lid, then press Read standard. When the reading is complete, remove Standard #2.
* Press Run samples.
* On the assay screen, select the sample volume and units
* Insert a sample tube into the sample chamber, close the lid, then press Read tube. When the reading is complete (~3 seconds), remove the sample tube.
* Repeat step last step until all samples have been read

## Visualize RNA on BioAnalyzer

#### Bioanalyzer Instructions

#### Prepping the Gel (use filtered gel within 4 weeks)

* Pipette 550uL gel matrix (red) into a spin filter
* Centrifuge at 1500g for 10 mins at room temp
* Aliquot 65uL filtered gel into 0.5mL RNase-free tubes

#### Prepping the Gel-Dye Mix

* Allow RNA dye (blue) to equilibrate to room temp for 30 mins (in the dark)
* Vortex RNA dye (blue) for 10 secs and spin down
* Add 1 uL RNA dye into a 65uL aliquot of filtered gel
* Vortex solution
* Spin gel-dye mix at 13,000g for 10 mins at room temp

#### Turn on Laptop/Bioanalyzer (use gel-dye mix within 1 day)

* Select correct icon
* Select correct program (RNA, DNA, protein)

#### Cleaning the Electrodes

* Slowly fill one of the wells of the electrode cleaner with 350uL RNase-free water
* Place the electrode cleaner in the bioanalyzer
* Close the lid and let sit for 5 mins
* Upon removing the electrode cleaner, keep lid open for at least 30 secs (to allow for evaporation)

#### Loading the Gel-Dye Mix

* Put a new chip on the chip priming station
* Pipette 9uL of gel-dye mix in the well marked G (white G on black background)
* Position the plunger at 1mL
* Close the chip priming station (MAKE SURE IT CLICKS)
* Press plunger down until it’s held by the clip
* Wait exactly 30 secs and release clip
* Wait exactly 5 secs and slowly pull plunger back to 1mL
* Open chip priming station and pipette 9uL gel-dye mix into wells marked G (black G on grey background)
* Discard remaining gel-dye mix

#### Loading the Conditioning Solution and Marker 

* Pipette 9uL of RNA conditioning solution (white) into well marked CS
* Pipette 5uL of RNA marker (green) in all 11 sample wells and the ladder well

#### Loading the Diluted Ladder and Samples

* Pipette 1uL ladder into the well marked ladder
* Pipette 1uL of sample into each of the 11 sample wells (pipette 1uL RNA marker (green) into each unused sample well)
* Vortex chip for 1 min at 2400 rpm
* Run the chip in the bioanalyzer within 5 mins (RNA assay)

**CLEAN THE ELECTRODES AGAIN AFTER EVERY RUN**

## Stranded mRNA-Seq Library Prep
### Using KAPA Stranded mRNA-Seq Kit using 1/2 rxn volumes
This should take 8-10 hours
Refer to manual during procedure (steps below are for notes and comments)

#### Reagents 

KAPA Stranded mRNA-Seq Kit (KAPA #KK8420). This kit includes all the enzymes and buffers required for cDNA library preparation from  isolation of poly(A)-tailed RNA. Kits include reagents for RNA fragmentation, 1st strand cDNA synthesis and 2nd strand synthesis/marking, and cDNA library preparation, including A-tailing, ligation and library amplification. 

Steps in Library construction:

* mRNA capture using magnetic oligo-dT beads
* Fragmentation using heat and magnesium
* 1st Strand cDNA Synthesis using random priming
* 2nd Strand cDNA Synthesis and marking, which converts the cDNA:RNA hybrid to double-stranded cDNA (dscDNA) and incorporates dUTP in the second cDNA strand
* A-tailing to add dAMP to the 3′-ends of the dscDNA library fragments
* Adapter ligation, where dsDNA adapters with 3′-dTMP overhangs are ligated to A-tailed library insert fragments
	* **NOTE** Here, we insert custom adapters.  See below.
* Library amplification to amplify library fragments carrying appropriate adapter sequences at both ends using high-fidelity, low-bias PCR; the strand marked with dUTP is not amplified.

#### Additional reagents needed:

Annealing buffer stock (10X):

| Component| Concentration|
|----------|--------------|
| Tris HCl, pH 8| 100 mM|
| NaCl|500 mM|
| EDTA| 10 mM|

10 mM Tris-HCL (pH 8.0 - 8.5)
#### Equipment 

* Magnetic stand and compatible tubes or striptubes 
* Thermocycler
* SPRI purification beads (KAPA Pure Beads or AmpureXP)

#### Custom Oligos needed to make adapters:

![alt text](/RNAadapter.png)


|Oligo Name| Sequence|
|----------|---------|
|Universal_SAI1_Adapter|AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCGTCGACT*T|
|Indexed_Adapter_SAI1_I5|P*AGTCGACGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG|
|Indexed_Adapter_SAI1_I8|P*AGTCGACGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG|
|Indexed_Adapter_SAI1_I9|P*AGTCGACGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG|
|Indexed_Adapter_SAI1_I11|P*AGTCGACGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG|


### Anneal RNA Adapters
Single-stranded oligos need to be annealed with their appropriate partner before ligation. 
 
* To create Adapter SAI1_I5, combine Universal_SAI1_Adapter with Indexed_Adapter_SAI1_I5 in a 1:1 ratio in working strength annealing buffer (final buffer concentration 1x) for a total annealed adapter concentration of 40uM (for example, if purchased oligos are resuspended to an initial concentration of 100uM, use 40uL Universal_SAI1_Adapter, 40ul Indexed_Adapter_SAI1_I5, 10ul 10x annealing buffer and 10ul nuclease-free water). Pair Universal_SAI1_Adapter with Indexed_Adapter_SAI1_I8, Indexed_Adapter_SAI1_I9, Indexed_Adapter_SAI1_I11 in the same fasion.
2.   In a thermocyler, incubate at 97.5°C for 2.5 minutes, and then cool at a rate of not greater than 3°C per minute until the solution reaches a temperature of 21°C. Hold at 4°C.
3.   Prepare final working strength concentrations of annealed adapters from this annealed stock (the appropriate working stock dilution for your experiment can be determined from the chart below). For convenience, it is possible to store the adapters at 4°C while in active use.  

#### Adapter concentration will vary depending on overall RNA yield, see table below:
|Quantity of starting material|Adapter stock concentration|Adapter concentration in ligation reaction|
|----|----|-----|
|100 – 250 ng|140 nM| 10 nM|
|251 – 500 ng|350 nM|25 nM|
|501 – 2000 ng|700 nM|50 nM|
|2001 – 4000 ng|1400 nM|100 nM|

For Puritz and Lotterhos 2017, we used 4000 ng starting RNA, but because of difficulties assessing and quantifying molluscan RNA, we chose to use a 700 nM working stock with a final reaction concentration of 50 nM.

**Procedure:**
* Prepare mastermixes for number of libraries (individual RNA extractions)
  * See tables in manual and tables below for guidelines (We are using 1/2 reactions):

|Component| Total Volume Needed for 4 RXNs (Includes 20% excess)|
|---------|--------------------|
|**1st Strand Synthesis Master Mix:**|--|
|1st Strand Synthesis Buffer|22 μl|
|KAPA Script|2 μl|
|**Total Master Mix Volume**| **24 μl**|
|**Final reaction composition:**|--|
|1st Strand Synthesis Master Mix| 5 μl|
|Fragmented, primed RNA|10 μl| 
|**Total Reaction Volume**| **15 μl**|

|Component| Total Volume Needed for 4 RXNs (Includes 10% excess)|
|---------|--------------------|
|**2nd Strand Synthesis and Marking Master Mix:**|--|
|2nd Strand Marking Buffer|62 μl|
|2nd Strand Synthesis Enzyme Mix|4 μl|
|**Total Master Mix Volume**| **66 μl**|
|**Final reaction composition:**|--|
|2nd Strand Synthesis and Marking Master Mix| 15 μl|
|Fragmented, primed RNA|15 μl| 
|**Total Reaction Volume**| **30 μl**|

|Component| Total Volume Needed for 4 RXNs (Includes 10% excess)|
|---------|--------------------|
|**A-Tailing Master Mix:**|--|
|Water|52.8 μl|
|10X KAPA A-Tailing Buffer|6.6 μl|
|KAPA A-Tailing Enzyme|6.6 μl|
|**Total Master Mix Volume**| **66 μl**|
|**Resuspend beads in a volume of:**|** 15 μl** |

|Component| Total Volume Needed for 4 RXNs (Includes 10% excess)|
|---------|--------------------|
|**Adapter Ligation Master Mix:**|--|
|Water|35.2 μl|
|5X KAPA Ligation Buffer| 30.8μl|
|KAPA T4 DNA Ligase|11 μl|
|**Total Master Mix Volume**| **77 μl**|
|**Final reaction composition:**|--|
|Beads with A-tailed DNA|15 μl|
|Adapter Ligation Master Mix|17.5 μl|
|Adapter (350 nM – 1400 nM, as appropriate)|2.5 μl|
|**Total Reaction Volume**| **35 μl**|

|Component| Total Volume Needed for 4 RXNs (Includes 10% excess)|
|---------|--------------------|
|**Library Ampli cation Master Mix:**|--|
|2X KAPA HiFi HotStart ReadyMix|55 μl|
|10X KAPA Library Amplication Primer Mix|11 μl|
|**Total Master Mix Volume**| **66 μl**|
|**Final reaction composition:**|--|
|Adapter-ligated library DNA|10 μl|
|Library Ampli cation Master Mix|15 μl|
|Balance of water (if required)|5 μl|
|**Total Reaction Volume**| **30 μl**|

### mRNA Capture

* Before mRNA capture beads can be used they must be washed with mRNA Bead Binding Buffer
	* Resuspend beads thoroughly by gentle pipetting or vortexing
	* For each library to be prepared, transfer 26.25 uL of the resuspended mRNA Capture beads into an appropriate tube
		* Up to 48 libraries (1,260 uL) can be washed in a single tube
	* Place the tube on a magnet holder and incubate at room temperature until solution is clear.
	* Discard supernatant and replace with an equal volume of mRNA Bead Binding Buffer.
	* Remove tube from magent and again resuspend the beads.
	* Place the tube on a magnet holder and incubate at room temperature until solution is clear.
	* Discard supernatant and replace with an equal volume of mRNA Bead Binding Buffer.
	* Remove tube from magent and again resuspend the beads.

* Combine the following for each RNA sample to be captured:

|Component|Volume|
|---------|------|
|RNA sample (in RNase-free water)| 25 μl|
|KAPA mRNA Capture Beads| 25 μl|
|**Total Volume**| **50 μl**|

* Mix thoroughly by gently pipetting up and down several times.
* Place the plate/tube in a thermal cycler and carry out the 1st mRNA capture program as follows:

|Step|Temp.|Duration|
|----|-----|--------|
|1st mRNA capture|65 °C|2 min|
|Cool|20 °C|5 min|


* Place the plate/tube containing the mixture of KAPA mRNA Capture Beads and RNA on a magnet and incubate at room temperature until the solution is clear. Remove and discard the supernatant.
* Remove the plate/tube from the magnet and resuspend thoroughly in 100 μl of KAPA mRNA Bead Wash Buffer by pipetting up and down several times.
* Place the plate/tube on the magnet and incubate at room temperature until the solution is clear. Remove and discard the supernatant.
* Resuspend the beads in 25 μl of RNase-free water.
* Place the plate/tube in a thermal cycler and carry out the 2nd mRNA capture program as follows:

|Step|Temp.|Duration|
|----|-----|--------|
|2nd mRNA capture|70 °C|2 min|
|Cool|20 °C|5 min|

* Add 25 μl of KAPA Bead Binding Buffer to the mixture of KAPA mRNA Capture Beads and RNA and mix thoroughly by gently pipetting up and down several times.
* Incubate the plate/tube at 20 °C for 5 min.
* Place the plate/tube on the magnet and incubate at room temperature until the solution is clear. Remove and discard the supernatant.
* Remove the beads from the magnet and resuspend in 100 μl of KAPA mRNA Bead Wash Buffer by pipetting up and down several times.
* Place the plate/tube on the magnet and incubate at room temperature until the solution is clear. Remove and discard the entire volume of supernatant.
#### mRNA Elution, Fragmentation and Priming
* Prepare the required volume of 1X Fragment, Prime and Elute Buffer as follows:

|Component|Volume per sample|
|---------|------|
|Water| 5.5 μl|
|Fragment, Prime and Elute Buffer (2X)| 5.5 μl|
|**Total Volume**| **11 μl**|

* Thoroughly resuspend the KAPA mRNA Capture Beads with captured mRNA in 11 μl of 1X Fragment, Prime and Elute Buffer.

---

### Safe Stopping Point
Resuspended beads with captured mRNA may be stored at 4 oC for up to 24 hours. Do not freeze the samples as this will damage the beads. When ready, proceed to step below.

---

* Place the plate/tubes in a thermal cycler and carry out the fragmentation and priming program as follows:

|Desired Fragment Size| Temp.| Duration|
|---------------------|------|---------|
|100 – 200 bp|94 °C|8 min|
|200 – 300 bp|94 °C|6 min|
|300 – 400 bp|85 °C|6 min|

##### For Puritz and Lotterhos (2017), we chose 94 °C for 7 mins to have fragments between 150-250 bp, approximately the same size distribution as planned for our DNA libraries.

* Immediately place the plate/tube on a magnet to capture the beads, and incubate until the liquid is clear. **Caution: To prevent hybridization of poly(A)- rich RNA to the capture beads, do not allow the sample to cool before placing on the magnet.**
* Carefully remove 10 μl of the supernatant containing the eluted, fragmented, and primed RNA into a separate plate or tube.
* Proceed immediately to **1st Strand Synthesis**.

### 1st Strand Synthesis
* On ice, assemble the 1st Strand Synthesis reaction as follows:

|Component|Volume|
|---------|------|
|Fragmented, primed RNA eluted from beads| 10 μl|
|1st Strand Synthesis Master Mix| 5 μl|
|**Total Volume**| **15 μl**|

* Keeping the plate/tube on ice, mix thoroughly by gently pipetting the reaction up and down several times.
* Incubate the plate/tube using the following protocol:

|Step|Temp.|Duration|
|----|-----|--------|
|Primer extension|25 °C|10 min|
|1st Strand synthesis|42 °C|15 min|
|Enzyme inactivation|70 °C|15 min|
|HOLD|4 °C|∞|

* Place the plate/tube on ice and proceed immediately to **2nd Strand Synthesis and Marking**.

### 2nd Strand Synthesis
* Assemble the 2nd Strand Synthesis and Marking reaction as follows:

|Component|Volume|
|---------|------|
|1st Strand cDNA| 15 μl|
|2nd Strand Synthesis and Marking Master Mix| 15 μl|
|**Total Volume**| **30 μl**|

* Mix thoroughly by gently pipetting the reaction up and down several times.
* Incubate the plate/tube using the following protocol:

|Step|Temp.|Duration|
|----|-----|--------|
|2nd Strand synthesis and marking|16 °C|60 min|
|HOLD|4 °C|∞|

* Place the plate/tube on ice and proceed immediately to **2nd Strand Synthesis and Marking Cleanup**.

### Cleanup

* Perform a 1.8X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|2nd Strand Synthesis reaction product| 30 μl|
|Agencourt® AMPure® XP reagent| 54 μl|
|**Total Volume**| **84 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 74 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. 
  * **Caution: over-drying the beads may result in dramatic yield loss.**
* Proceed immediately to **A-Tailing** immediately, or follow the Safe Stopping Point instructions below.

---

### SAFE STOPPING POINT
* Resuspend the beads in 15 μl 1X A-Tailing Buffer (Table 5B), cover the reaction and store at 4 oC for up to 24 hours. Do not freeze the samples as this will damage the AMPure® XP® beads. When ready, proceed to **A-Tailing after Safe Stopping Point**.

---

### A-Tailing
* A-Tailing is performed either directly after the 2nd Strand Synthesis and Marking Cleanup, or after the Safe Stopping Point, where beads were resuspended in 1X A-Tailing Buffer and stored at 4 °C for up to 24 hours. 

#### A-Tailing immediately
* Assemble the A-Tailing reaction as follows:

|Component|Volume|
|---------|------|
|Beads with dscDNA| --|
|A-Tailing Master Mix| 15 μl|
|**Total Volume**| **15 μl**|

* Mix thoroughly by pipetting up and down several times.
* Incubate the plate/tube using the following protocol:

|Step|Temp.|Duration|
|----|-----|--------|
|A-Tailing|30 °C|30 min|
|Enzyme inactivation|60 °C|30 min|
|HOLD|4 °C|∞|

* Proceed immediately to **Adapter Ligation**.

#### A-Tailing after safe stopping point
* To resume library preparation, combine the following reagents to perform A-Tailing:

|Component|Volume|
|---------|------|
|Beads with dscDNA (in 1X A-Tailing Buffer)| 7.5 μl |
|A-Tailing Master Mix after Safe Stopping Point| 7.5 μl|
|**Total Volume**| **15 μl**|

* Mix thoroughly by pipetting up and down several times.
* Incubate the plate/tube using the following protocol:

|Step|Temp.|Duration|
|----|-----|--------|
|A-Tailing|30 °C|30 min|
|Enzyme inactivation|60 °C|30 min|
|HOLD|4 °C|∞|

* Proceed immediately to **Adapter Ligation**.

### Adapter Ligation

#### Adapter concentration will vary depending on overall RNA yield, see table below:

|Quantity of starting material|Adapter stock concentration|Adapter concentration in ligation reaction|
|----|----|-----|
|100 – 250 ng|140 nM| 10 nM|
|251 – 500 ng|350 nM|25 nM|
|501 – 2000 ng|700 nM|50 nM|
|2001 – 4000 ng|1400 nM|100 nM|

For Puritz and Lotterhos 2017, we used 4000 ng starting RNA, but because of difficulties assessing and quantifying molluscan RNA, we chose to use a 700 nM working stock with a final reaction concentration of 50 nM.

#### This will be where we insert the custom adapters that are barcoded with RE sites

* Set up the adapter ligation reactions as follows:

|Component|Volume|
|---------|------|
|Beads with A-tailed DNA| 15 μl |
|Adapter Ligation Master Mix| 17.5 μl |
|**Adapters***| 2.5 μl|
|**Total Volume**| **35 μl**|


* Mix thoroughly by pipetting up and down several times to resuspend the beads.
* Incubate the plate/tube at 20 °C for 30 min.
* Proceed immediately to **1st Post-Ligation Cleanup**.

### Post-Ligation Cleanup

* Perform a 1X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|Beads with adapter-ligated DNA| 35 μl|
|Agencourt® AMPure® XP reagent| 35 μl|
|**Total Volume**| **70 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 65 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the beads in 25 μl of 10 mM Tris-HCl (pH 8.0).
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads

---

### Safe Stopping Point
The solution with resuspended beads can be stored at 4 °C for up to 24 hours. Do not freeze the beads, as this can result in dramatic loss of DNA. When ready, proceed to **2nd Post-Ligation Cleanup**.

---

### 2nd Post-Ligation Cleanup

* Perform a 1X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|Beads with puri ed, adapter-ligated DNA| 25 μl|
|Agencourt® AMPure® XP reagent| 25 μl|
|**Total Volume**| **50 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 45 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the beads in 11.25 μl of 10 mM Tris-HCl (pH 8.0).
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 20 μl of the clear supernatant to a new plate/tube and proceed to *Library Amplication*.
---

### SAFE STOPPING POINT
The purified, adapter-ligated library DNA may be stored at 4 °C for up to 1 week, or frozen at -20 °C for up to 1 month. When ready, proceed to **Library Amplification**.

---
### Library Amplificiation

* Assemble each library ampli cation reaction as follows:

|Component|Volume|
|---------|------|
|Purified, adapter-ligated DNA| 10 μl|
|Library Amplification Master Mix| 15 μl|
|**Total Volume**| **25 μl**|

* Mix well by pipetting up and down several times
* Amplify the library using the following thermal cycling protocol:

|Step|Temp|Duration|Cycles|
|----|----|--------|------|
|Initial denaturation|98 °C|45 sec|1|
|Denaturation|98 °C|15 sec|12|
|Annealing*|60 °C|30 sec|12|
|Extension|72 °C|30 sec|12|
|Final Extension|72 °C|5 min|1|
|Hold|10 °C | ∞|1|


* Place the plate/tube on ice and proceed to **Library Amplification Cleanup**

###Library Amplification Cleanup

* Perform a 1X SPRI® cleanup by combining the following

|Component|Volume|
|---------|------|
|Amplified library DNA| 25 μl|
|Agencourt® AMPure® XP reagent| 25 μl|
|**Total Volume**| **50 μl**|

* Mix thoroughly by pipetting up and down several times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 45 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the dried beads in 22 μl of 10 mM Tris-HCl (pH 8.0).
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
Transfer 20 μl of the clear supernatant to a new plate/tube.

## Quant libraries
**Procedure (Standard HS DNA protocol)**
* Set up the required number of 0.5-mL tubes for standards and samples. The Qubit® RNA HS Assay requires 2 standards.
* Label the tube lids.
* Prepare the Qubit® working solution by diluting the Qubit® DNA HS Reagent 1:200 in Qubit® DNA HS Buffer. Use a clean plastic tube each time you prepare Qubit® working solution. **Do not mix the working solution in a glass container.**
* Add 190 μL of Qubit® working solution to each of the tubes used for standards.
* Add 10 μL of each Qubit® standard to the appropriate tube, then mix by vortexing 2–3 seconds. Be careful not to create bubbles.
* Add Qubit® working solution to individual assay tubes so that the final volume in each tube after adding sample is 200 μL.
* Add each sample to the assay tubes containing the correct volume of Qubit® working solution, then mix by vortexing 2–3 seconds. The final volume in each tube should be 200 μL.
* Allow all tubes to incubate at room temperature for 2 minutes.
* On the Home screen of the Qubit® 3.0 Fluorometer, press DNA, then select DNA: High Sensitivity as the assay type. The “Read standards” screen is displayed. Press Read Standards to proceed.
* Insert the tube containing Standard #1 into the sample chamber, close the lid, then press Read standard. When the reading is complete (~3 seconds), remove Standard #1.
* Insert the tube containing Standard #2 into the sample chamber, close the lid, then press Read standard. When the reading is complete, remove Standard #2.
* Press Run samples.
* On the assay screen, select the sample volume and units
* Insert a sample tube into the sample chamber, close the lid, then press Read tube. When the reading is complete (~3 seconds), remove the sample tube.
* Repeat step last step until all samples have been read

---

### Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---


### DSN Normalization
#### DSN needs to be properly dilued and should be tested for activity levels before proceeding

#### The protocol below was taken from Illumina's recommendations [LINK](https://support.illumina.com/content/dam/illumina-support/documents/myillumina/7836bd3e-3358-4834-b2f7-80f80acb4e3f/dsn_normalization_sampleprep_application_note_15014673_c.pdf)
#### Reagents

| Reagent| Supplier|
|----------|--------------|
|1 M HEPES buffer solution|Invitrogen, part # 15630‐080 |
|5 M NaCl solution|Ambion, part # AM9760G|
|KAPA HiFi HotStart PCR kit with dNTPs|Kapa, part #KK2502|
|Strip tubes|General lab supplier|
|DSN Kit|Evrogen, part # EA001 Sigma Aldrich, part # E7023|
|Ethanol 200 proof (absolute) for molecular biology (500 ml)|AB, part # 4333764F|
|PCR Primer PE 1.0|Included in Kapa stranded mRNA kit|
|PCR Primer PE 2.0|Included in Kapa stranded mRNA kit|
|SPRI beads|Agencourt AMPure, part # 29152; KAPA Pure Beads, part #KK8000|
|Nuclease-free water|General lab supplier|

#### Equipment
* Thermocycler
* Magentic stand compatible with strip tubes

#### Procedure

* First pool individual RNA libraries in equal quantities to create a single pool of 500 ng.
	* For example pool 125 ng each of four individual libraries.

* Create a 4X hybridization solution

|Component|Volume|
|---------|------|
|1 M HEPES buffer solution| 200 𝜇l|
|5 M NaCl solution| 400 𝜇l|
|Nuclease‐free water| 400 𝜇L|
|**Total Volume**|**1000 𝜇**L|

* Use two thermocyclers and set one to hold at 68°C
* Prepare the following reaction mix in a separate, sterile, nuclease‐free 200 μl PCR tube on ice for each sample to be normalized.

|Component|Volume|
|---------|------|
|Sample library (500 ng)| 13.5 𝜇l|
|4X Hybridization buffer| 4.5 𝜇l|
|**Total Volume Per Sample**|**18 𝜇**L|

* Gently pipette the entire volume up and down 10 times, then centrifuge briefly to mix.
* Transfer the entire volume of reaction mix directly to the bottom of a new, sterile, nuclease‐free 200 μl PCR tube, using a pipette. Do not let the sample contact the side of the tube during the process.
* Incubate the reaction mix tube on the thermal cycler using the following PCR cycling conditions:

|Step|Temp|Duration|
|----|----|--------|
|Initial denaturation|98 °C|2 min|
|Treatment|68 °C|5 hours|

* **Caution**- Following incubation, keep the thermal cycler lid closed and the temperature held at 68°C. Do not remove the reaction mix tube from thermal cycler prior to and during DSN treatment.
* Dilute the 10X DSN Master buffer supplied in the DSN kit to 2X with nuclease‐ free water
* Pre‐heat the 2X DSN buffer on the pre‐heated heat block at 68°C.
  * **Note:** Do not remove the 2X DSN buffer from the heat block during DSN treatment.
* Quickly add 20 μl of pre‐heated 2X DSN buffer to the first reaction mix tube.
* With the reaction mix tube remaining within the thermal cycler, gently pipette the entire volume up and down 10 times to mix thoroughly using a pipette set to 40 μl.
  * **Note**:Pipette the solution directly to the bottom of the PCR tube and do not let the sample contact the side of the tube during the process.
  * **Note**: It is important to keep the thermal cycler closed, except for the time necessary to add the 2X DSN buffer and mix. When preparing more than one reaction mix tube, keep the thermal cycler lid closed when extracting the 2X DSN buffer from its tube, then open the thermal cycler lid only for the time necessary to add and mix the 2X DSN buffer.

* Repeat steps 2 and 3 for each reaction mix tube.
* Incubate the reaction mix tubes on the thermal cycler at 68°C for 10 minutes.
* Quickly add 2 μl of DSN enzyme to the first reaction mix tube using a 2 μl pipette.
* With the reaction mix tube remaining within the thermal cycler, gently pipette the entire volume up and down 10 times to mix thoroughly using a pipette set to 40 μl.
  * **Note**:Pipette the solution directly to the bottom of the PCR tube and do not let the sample contact the side of the tube during the process.
* Repeat steps 6 and 7 for each reaction mix tube.
* Incubate the reaction mix tubes on the thermal cycler at 68°C for 25 minutes.
* Add 40 μl of 2X DSN stop solution to each reaction mix tube. Gently pipette the entire volume up and down to mix thoroughly, then place the tubes on ice.

---

### Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---
### SPRI Cleanup

* Perform a 1X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|DSN Treated Library| 80 μl|
|Agencourt® AMPure® XP reagent| 160 μl|
|**Total Volume**| **240 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 235 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the beads in 25 μl of 10 mM Tris-HCl (pH 8.0).
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 24 μl of the clear supernatant to a new plate/tube and proceed to next step.

### PCR Enrichment

|Component|Volume|
|---------|------|
|DSN Treated Library | 30 μl|
|2X KAPA HiFi HotStart ReadyMix| 25 μl|
|10X KAPA Library Amplification Primer Mix| 5 μl|
|**Total Volume per sample**| **50 μl**|

* *For now, let's use the same master mix from the Kappa Kit*
* Mix well by pipetting up and down several times
* Amplify the library using the following thermal cycling protocol:

|Step|Temp|Duration|Cycles|
|----|----|--------|------|
|Initial denaturation|98 °C|45 sec|1|
|Denaturation|98 °C|15 sec|12|
|Annealing*|60 °C|30 sec|12|
|Extension|72 °C|30 sec|12|
|Final Extension|72 °C|5 min|1|
|Hold|10 °C | ∞|1|

### SPRI Cleanup

* Perform a 1.6X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|Enriched DSN Library| 50 μl|
|Agencourt® AMPure® XP reagent| 80 μl|
|**Total Volume**| **130 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 115 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the beads in 22 μl of 10 mM Tris-HCl (pH 8.0).
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 22 μl of the clear supernatant to a new plate/tube and proceed to next step.

## Quant libraries
**Procedure (Standard HS DNA protocol)**
* Set up the required number of 0.5-mL tubes for standards and samples. The Qubit® RNA HS Assay requires 2 standards.
* Label the tube lids.
* Prepare the Qubit® working solution by diluting the Qubit® DNA HS Reagent 1:200 in Qubit® DNA HS Buffer. Use a clean plastic tube each time you prepare Qubit® working solution. **Do not mix the working solution in a glass container.**
* Add 190 μL of Qubit® working solution to each of the tubes used for standards.
* Add 10 μL of each Qubit® standard to the appropriate tube, then mix by vortexing 2–3 seconds. Be careful not to create bubbles.
* Add Qubit® working solution to individual assay tubes so that the final volume in each tube after adding sample is 200 μL.
* Add each sample to the assay tubes containing the correct volume of Qubit® working solution, then mix by vortexing 2–3 seconds. The final volume in each tube should be 200 μL.
* Allow all tubes to incubate at room temperature for 2 minutes.
* On the Home screen of the Qubit® 3.0 Fluorometer, press DNA, then select DNA: High Sensitivity as the assay type. The “Read standards” screen is displayed. Press Read Standards to proceed.
* Insert the tube containing Standard #1 into the sample chamber, close the lid, then press Read standard. When the reading is complete (~3 seconds), remove Standard #1.
* Insert the tube containing Standard #2 into the sample chamber, close the lid, then press Read standard. When the reading is complete, remove Standard #2.
* Press Run samples.
* On the assay screen, select the sample volume and units
* Insert a sample tube into the sample chamber, close the lid, then press Read tube. When the reading is complete (~3 seconds), remove the sample tube.
* Repeat step last step until all samples have been read

---

## Split finished cDNA library for each sample into two vials
* One tube for sequencing
* One tube for probe synthesis

---

### Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---

## Probe Synthesis

Review quantifications for probes.  Ideally, there should be about 500 ng of probes per capture.  If not enough of probes is obtained, the PCR product can be re-amplified.

### Remove adapters from cDNA

#### Materials needed
| reagent                                | producer         | serial #      |            
|----------------------------------------|------------------|-----------|
|Mung Bean Nuclease| NEB| M0250S|
|SalI-HF| NEB| R3138T|
|Agencourt AMPure XP  |Beckman Coulter   | A63881|

* Setup a restriction digest using 1 μg of DSN library

|Component|Volume|
|---------|------|
| 10X Cutsmart Buffer| 4 μl|
|SalI-HF Enzyme (100 units)| 1 μl |
| Molecular Grade H20| 22.75 μl|
| DSN Enriched Library| 12.25 μl|
|**Total Volume**| **40 μl**|

* Incubate reactions in thermocycler at 37°C for at least 8 hours, prefereably 12-16 hours.

* **Using the same tubes from the previous step** combine the following:

|Component|Volume|
|---------|------|
|Digested Library| 40 μl|
|10X Mung Bean Nuclease buffer| 4.5 μl|
|Mung Bean Nuclease (10 units per μl)| 0.5 μl|
|**Total Volume**| **45 μl**|

* Incubate at 30°C for 30 minutes

* Perform a 1.8X SPRI cleanup by combining the following:

#### This step may be possbile to skip and proceed directly to the 1.5X SPRI Cleanup.

|Component|Volume|
|---------|------|
|MBN reaction| 45 μl|
|Agencourt® AMPure® XP reagent| 81 μl|
|**Total Volume**| **126 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 115 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the beads in 22 μl of 10 mM Tris-HCl (pH 8.0). Volume needed depends on the number of captures. Calculate 10 μl per capture plus an aliquot for checking the probes concentration using Qubit.
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 20 μl of the clear supernatant to a new plate/tube and proceed to next step.

---

### Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---

* Perform a 1.5X SPRI cleanup by combining the following:

|Component|Volume|
|---------|------|
|MBN reaction| 22 μl|
|Agencourt® AMPure® XP reagent|33 μl|
|**Total Volume**| **55 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 115 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the beads in 31 μl of 10 mM Tris-HCl (pH 8.0). Volume needed depends on the number of captures. Calculate 10 μl per capture plus an aliquot for checking the probes concentration using Qubit.
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 30 μl of the clear supernatant to a new plate/tube and proceed to next step.

---

### Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---


##Biotin Labeling

###Materials needed

| reagent                                | producer         | serial #      |            
|----------------------------------------|------------------|-----------|
|DecaLabel™ Biotin DNA Labeling Kit |Thermo Scientific | FERK0651   |                     

###Procedure

* Add the following components into 1.5 ml microcentrifuge tube:

|Component|Volume|
|---------|------|
|RE treated DSN Library| 20 μl|
|Decanucleotide in 5X Reaction Buffer| 10 μl|
| Water, nuclease-free| 14 μl|
|**Total Volume**| **44 μl**|

* Vortex the tube and spin down in a microcentrifuge for 3-5 s
* Incubate the tube in a boiling water bath for 5-10 min and cool it on ice. Spin down quickly.
* Add the following components in the same tube:

|Component|Volume|
|---------|------|
|Biotin Labeling Mix| 5 μl|
|Klenow fragment, exo– (5 u)| 1 μl|
|**Total Volume**| **50 μl**|

* Shake the tube and spin down in a microcentrifuge for 3-5 s. 
* Incubate for 12-20 hours at 37°C. 

####Optional:Control reaction
* Add the following components into 1.5 ml microcentrifuge tube:

|Component|Volume|
|---------|------|
|Control Template, 10 ng/μl | 25 μl|
|Decanucleotide in 5X Reaction Buffer| 10 μl|
| Water, nuclease-free| 9 μl|
|**Total Volume**| **44 μl**|

* Vortex the tube and spin down in a microcentrifuge for 3-5 s
* Incubate the tube in a boiling water bath for 5-10 min and cool it on ice. Spin down quickly.
* Add the following components in the same tube:

|Component|Volume|
|---------|------|
|Biotin Labeling Mix| 5 μl|
|Klenow fragment, exo– (5 u)| 1 μl|
|**Total Volume**| **50 μl**|


* Perform a 1.5X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|Biotin reaction| 50 μl|
|Agencourt® AMPure® XP reagent| 75 μl|
|**Total Volume**| **125 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 115 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the beads in 20 μl of 10 mM Tris-HCl (pH 8.0). Volume needed depends on the number of captures. Calculate 10 μl per capture plus an aliquot for checking the probes concentration using Qubit.
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 19 μl of the clear supernatant to a new plate/tube and proceed to next step.

---

###Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---


##Preparation of whole genome libraries using KAPA HyperPlus Kit

Refer to [manual](https://www.kapabiosystems.com/document/kapa-hyperplus-library-preparation-kit-tds/?dl=1) during procedure (steps below are for notes and comments). 

#### This assumes that genomic DNA is already extracted and sheared.


Oligos needed:

|Name| 5' to 3' Sequence|
|----------|---------------------------------------------|
|DNA_P1.1.1|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCATGG*T|
|DNA_P1.1.2|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTAACCAG*T|
|DNA_P1.1.3|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTCGATCG*T|
|DNA_P1.1.4|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTTCGATG*T|
|DNA_P1.1.5|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTTGCATG*T|
|DNA_P1.1.6|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTCAACCG*T|
|DNA_P1.1.7|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTGGTTGG*T|
|DNA_P1.1.8|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTAAGGAG*T|
|DNA_P1.1.9|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTAGCTAG*T|
|DNA_P1.1.10 |	ACACTCTTTCCCTACACGACGCTCTTCCGATCTACACAG*T|
|DNA_P1.1.11|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTAATTAG*T|
|DNA_P1.1.12|	ACACTCTTTCCCTACACGACGCTCTTCCGATCTACGGTG*T|
|DNA_P1.2.1|	PC*CGTACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.2|	PC*TTGGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.3|	PC*GCTAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.4|	PC*AGCTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.5|	PC*ACGTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.6|	PC*GTTGGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.7|	PC*CCAACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.8|	PC*TTCCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.9|	PC*TCGATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.10|	PC*TGTGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.11|	PC*TTAATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P1.2.12|	PC*TGCCAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
|DNA_P2.1|	P*GATCGGAAGAGCGAGAACAA|
|DNA_P2.2|	GTGACTGGAGTTCACACGTGTGCTCTTCCGATC*T|

### Anneal Adapters
Single-stranded oligos need to be annealed with their appropriate partner before ligation. 
 
* To create Adapter P1, combine each oligo 1.1 with its complementary oligo 1.2 in a 1:1 ratio in working strength annealing buffer (final buffer concentration 1x) for a total annealed adapter concentration of 40uM (for example, if purchased oligos are resuspended to an initial concentration of 100uM, use 40ul oligo 1.1, 40ul oligo 1.2, 10ul 10x annealing buffer and 10ul nuclease-free water). Do the same for oligos 2.1 and 2.2 to create the common adapter P2. 
2.   In a thermocyler, incubate at 97.5°C for 2.5 minutes, and then cool at a rate of not greater than 3°C per minute until the solution reaches a temperature of 21°C. Hold at 4°C.
3.   Prepare final working strength concentrations of annealed adapters from this annealed stock (the appropriate working stock dilution for your experiment can be determined from our ligation molarity calculator). For convenience, it is possible to store the adapters at 4°C while in active use.


### End repair
* Adjust sample volume for fragmented DNA samples to 25 μl.
* Add the following to each sample:

|Component|Volume|
|---------|------|
|End Repair & A-Tailing Buffer* | 3.5 μl|
|End Repair & A-Tailing Enzyme Mix* | 1.5 μl|
|Fragmented, double-stranded DNA| 25 μl |
|**Total Volume**| **30 μl**|

  * The buffer and enzyme mix should preferably be pre-mixed and added in a single pipetting step. 
    * Premixes are stable for ≤24 hrs at room temperature, for ≤3 days at 4°C, and for ≤4 weeks at -20°C
  * Vortex gently and spin down brie y. Return the reaction plate/tube(s) to ice. Proceed immediately to the next step.
  * Incubate in a thermocycler programmed as outlined below. A heated lid is required for this step. If possible, set the temperature of the heated lid to ~85°C (instead of the usual 105°C).

|Step|Temp|Time|
|----|----|--------|
|End repair and A-tailing 1|20 °C|30 min|
|End repair and A-tailing 2|65 °C|30 min|
|Hold|10 °C | ∞|

* Notes
  * A heated lid is required for this incubation. If possible, set the temperature of the lid at 85°C, instead of the usual ~105°C.
  * If proceeding to the adapter ligation reaction setup without any delay, the reaction may be cooled to 20°C instead of 4°C.

### Adapter ligation
* Dilute adapter stocks to the appropriate concentration, as outlined below:

|Fragmented DNA| Adapter stock concentration | Adapter:insert molar ratio|
|--------------|-----------------------------|---------------------------|
| 1 μg|15 μM|10:1|
|500 ng|15 μM|20:1|
|250 ng|15 μM|40:1|
|100 ng|15 μM|100:1|
|50 ng|15 μM|200:1|
|25 ng|15 μM|200:1|
|10 ng|15 μM|200:1|
|5 ng|15 μM|200:1|
|2.5 ng|15 μM|200:1|
|1 ng|15 μM|200:1|

* In the same plate/tube(s) in which end repair and A-tailing was performed, assemble each adapter ligation reaction as follows:

|Component|Volume|
|---------|------|
|End repair and A-tailing reaction product| 30 μl|
|Adapter stock (concentration as required) | 2.5 μl|
|PCR-grade water*| 2.5 μl |
|Ligation Buffer*| 30 μl |
|DNA ligase*| 5 μl |
|**Total Volume**| **55 μl**|

* Notes
  * The water, buffer and ligase enzyme should preferably be premixed and added in a single pipetting step. Premixes are stable for ≤24 hrs at room temperature, for ≤3 days at 4°C, and for ≤4 weeks at -20°C.
  
* Mix thoroughly and centrifuge briefly.
* Incubate at 20°C for 15 min.
  * Note: to achieve higher conversion rates and library yields, particularly for low-input samples, consider increasing the ligation time—to a maximum of 4 hrs at 20°C, or overnight at 4°C. Please note that longer ligation times may lead to increased levels of adapter-dimer. Adapter concentrations may have to be optimized if ligation times are extended signi cantly.
* Proceed immediately to the next step.

### Post-ligation Cleanup
* In the same plate/tube(s), perform a 0.8X bead- based cleanup by combining the following:

|Component|Volume|
|---------|------|
|Adapter ligation reaction product| 55 μl|
|KAPA Pure Beads | 44 μl|
|**Total Volume**| **99 μl**| 

* Mix thoroughly by vortexing and/or pipetting up and down multiple times.
* Incubate the plate/tube(s) at room temperature for 5 – 15 min to bind DNA to the beads.
* Place the plate/tube(s) on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard the supernatant.
* Keeping the plate/tube(s) on the magnet, add 200 μL of 80% ethanol.
* Incubate the plate/tube(s) on the magnet at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube(s) on the magnet, add 200 μL of 80% ethanol.
* Incubate the plate/tube(s) on the magnet at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature for 3 – 5 min, or until all of the ethanol has evaporated. *Caution: over-drying the beads may result in reduced yield.*
* Remove the plate/tube(s) from the magnet.
* Thoroughly resuspend the beads in in 12 μL of elution buffer (10 mM Tris-HCl, pH 8.0 – 8.5)
* Incubate the plate/tube(s) at room temperature for 2 min to elute DNA off the beads.
* Place the plate/tube(s) on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 11 μL of supernatant to a new plate/tube(s):

### Quant samples
**Procedure (Standard HS DNA protocol)**
* Set up the required number of 0.5-mL tubes for standards and samples. The Qubit® RNA HS Assay requires 2 standards.
* Label the tube lids.
* Prepare the Qubit® working solution by diluting the Qubit® DNA HS Reagent 1:200 in Qubit® DNA HS Buffer. Use a clean plastic tube each time you prepare Qubit® working solution. **Do not mix the working solution in a glass container.**
* Add 190 μL of Qubit® working solution to each of the tubes used for standards.
* Add 10 μL of each Qubit® standard to the appropriate tube, then mix by vortexing 2–3 seconds. Be careful not to create bubbles.
* Add Qubit® working solution to individual assay tubes so that the final volume in each tube after adding sample is 200 μL.
* Add each sample to the assay tubes containing the correct volume of Qubit® working solution, then mix by vortexing 2–3 seconds. The final volume in each tube should be 200 μL.
* Allow all tubes to incubate at room temperature for 2 minutes.
* On the Home screen of the Qubit® 3.0 Fluorometer, press DNA, then select DNA: High Sensitivity as the assay type. The “Read standards” screen is displayed. Press Read Standards to proceed.
* Insert the tube containing Standard #1 into the sample chamber, close the lid, then press Read standard. When the reading is complete (~3 seconds), remove Standard #1.
* Insert the tube containing Standard #2 into the sample chamber, close the lid, then press Read standard. When the reading is complete, remove Standard #2.
* Press Run samples.
* On the assay screen, select the sample volume and units
* Insert a sample tube into the sample chamber, close the lid, then press Read tube. When the reading is complete (~3 seconds), remove the sample tube.
* Repeat step last step until all samples have been read

###Library Amplification

* Assemble each library ampli cation reaction as follows:

|Component|Volume|
|---------|------|
|KAPA HiFi HotStart ReadyMix (2X) | 12.5 μl|
|Index PCR Primers F+R  | 2.5 μl|
|Adapter-ligated library| 10.0 μl|
|**Total Volume**| **25 μl**| 

* Calculate number of cycles needed based on previous quants

|Amount of adapter- ligated DNA in ampli cation reaction| Number of cycles required to generate 1 μg of library DNA| 
|--------------|-------------------------------------------|
|500 ng|1-2|
|100 ng|3-4|
|50 ng|5-6|
|25 ng|7-8|
|10 ng|8-9|
|5 ng|11-12|
|1 ng|12-13|

* Mix thoroughly and centrifuge briefly.
* Amplify using the following cycling protocol:

|Step|Temp|Duration|Cycles|
|----|----|--------|------|
|Initial denaturation|98 °C|45 sec|1|
|Denaturation|98 °C|15 sec|X|
|Annealing*|60 °C|30 sec|X|
|Extension|72 °C|30 sec|X|
|Final Extension|72 °C|1 min|1|
|Hold|4 °C | ∞|1|

* Proceed immediately to the next step

### Post-amplification Cleanup

* In the library amplification plate/tube(s) perform a 1X bead-based cleanup by combining the following:

|Component|Volume|
|---------|------|
|Adapter ligation reaction product| 25 μl|
|KAPA Pure Beads | 25 μl|
|**Total Volume**| **50 μl**| 

* Mix thoroughly by vortexing and/or pipetting up and down multiple times
* Incubate the plate/tube(s) at room temperature for 5 – 15 min to bind DNA to the beads
* Place the plate/tube(s) on a magnet to capture the beads. Incubate until the liquid is clear
* Carefully remove and discard the supernatant
* Keeping the plate/tube(s) on the magnet, add 200 μL of 80% ethanol.
* Incubate the plate/tube(s) on the magnet at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube(s) on the magnet, add 200 μL of 80% ethanol.
* Incubate the plate/tube(s) on the magnet at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature for 3 – 5 min, or until all of the ethanol has evaporated. Caution: over-drying the beads may result in reduced yield.
* Remove the plate/tube(s) from the magnet.
* Resuspend in 15 μl of 10 mM Tris or water

---

###Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---

## Hybridization and Capture

####Materials needed
| reagent                                | producer         | serial #      |            
|----------------------------------------|------------------|-----------|
|Denhardt’s solution  (50x)              |Life Technologies | 750018        |           
|Dynabeads® M-280 Streptavidin           |Life Technologies | 11205D, M-270 |           
|SSC Buffer Concentrate (20x)            |Fisher Scientific | 5075059      |              
|SDS Micropellets           |Fisher Scientific | BP8200100      |             
|Cot-1 DNA (1 mg/ml)                     |ThermoFischer     | 15279011|                   
|Agencourt AMPure XP                     |Beckman Coulter   | A63881|


|name       |5' to 3' Sequence                                                |
|-----------|------------------------------------------------------------------|
| BO1.P5.F  | AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT       |
| BO2.P5.R  | AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT       |
| BO3.P7.F  | CAAGCAGAAGACGGCATACGAGATIIIIIIGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT |
| BO4.P7.R  | AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIATCTCGTATGCCGTCTTCTGCTTG |


Solutions to prepare:

- 10 mM Tris-HCl pH 8.5 or PCR-grade water
- EDTA 500 mM
- SDS 10%
- TEN (10 mM Tris-HCl pH 7.5, 1 mM EDTA, 1M NaCl)
- 1x SSC / 0.1% SDS
- 0.5x SSC / 0.1% SDS
- 0.1x SSC / 0.1% SDS

Protocol based on previously described methods [hyRAD](https://github.com/chiasto/hyRAD/blob/master/wetlab.md#4-hybridization-capture-and-library-re-amplification) and [general capture](http://openwetware.org/wiki/Hyb_Seq_Prep)

Test the tubes and caps used for evaporation before starting by filling with known amount of water and keeping at 65°C for 48 hours.

Remember to perform one capture per pool of the libraries amplified with the same Illumina indexed primer.

###Hybridization

* Prepare the hybridization mix. Probes and blocking oligos are used in excess. Let's aim for 500 ng of probes and library

| Component                                | Volume  |
| -------------------------------------- | ----------- |
| water                                  | 3.5 μl        |
| SSC (20x)                              | 12.0  μl      |
| EDTA (500 mM)                          | 0.4 μl        |
| SDS (10%)                              | 0.4 μl        |
| Denhardt’s solution (50x)              | 1.6 μl        |
| Cot-1 DNA (1 mg/ml)                    | 0.5 μl        |
| BO.1 blocking oligo (200 μM)           | 0.4 μl        |
| BO.2 blocking oligo (200 μM)           | 0.4 μl        |
| BO.3 blocking oligo (200 μM)           | 0.4 μl        |
| BO.4 blocking oligo (200 μM)           | 0.4 μl        |
| prepared Illumina library (500 ng)| 10.0 μl        |
| probes (500 ng)               | 10.0 μl         |


* Incubate at 95°C for 10 minutes, then at 65°C for 48 hours. Mix from time to time. The best is to use rotor in the hybridization oven, but it is also possible to use a standard PCR machine.

### Preparation of Dynabeads

* Resuspend well Dynabeads M-280 (10 mg/ml).
2. Dispense 10 μl of beads in a PCR tube.
3. Wash:
    - magnetize, remove and discard supernatant,
    - resuspend in 200 μl of TEN.

4. Perform previous step 3 times in total.
5. Store in 200 μl of TEN at room temperature until use.
  * If more captures are expected, increase the initial amount of beads accordingly, transfer the final resupension into an eppendorf tube and add the appropriate volume of TEN (10 μl of beads should be resuspended in 200 μl of TEN).

### Washes

* Add 40 μl of the hybridization mixture to 200ul of Dynabeads
* Gently mix with pippette or inverting tube
* Incubate 30 min at room temperature.	


* Place on the magnet
* Remove supernatant and retain in case of DNA loss.
* Resuspend beads in 200 μl of 65°C 1x SSC / 0.1% SDS. 
* Mix well and incubate for 15 min, 65°C.	


* Place on the magnet
* Remove supernatant and retain in case of DNA loss.
* Replace with 200 μl of 65°C 1x SSC / 0.1% SDS. 
* Mix well and incubate for 10 min, 65°C.	


* Place on the magnet
* Remove supernatant and retain in case of DNA loss.
* Replace with 200 μl of 0.5x SSC / 0.1% SDS.
* Mix well, incubate for 10 min, 65°C.	
	
* Place on the magnet
* Remove supernatant and retain in case of DNA loss.
* Replace with 200 μl of 0.1x SSC / 0.1% SDS.
* Mix well, incubate for 10 min, 65°C.	
	
* Place on the magnet
* Remove supernatant and retain in case of DNA loss.
* Replace with 30 μl of 80°C H20 
* Mix well, incubate for for 10 min, 80°C.
* Place on magnet
* Remove and **SAVE** supernatant (this contains the hybridization-enriched products)
* Discard the beads.						

### Library re-amplification

* Assemble each library ampli cation reaction as follows:

|Component|Volume|
|---------|------|
|KAPA HiFi HotStart ReadyMix (2X) | 12.5 μl|
|Index PCR Primers F+R  | 2.5 μl|
|Enriched Library| 10.0 μl|
|**Total Volume**| **25 μl**| 

**NOTE:*** It's important to use the same INDEX primer as the original library prep!

* Mix thoroughly and centrifuge briefly.
* Amplify using the following cycling protocol:

|Step|Temp|Duration|Cycles|
|----|----|--------|------|
|Initial denaturation|98 °C|45 sec|1|
|Denaturation|98 °C|15 sec|12|
|Annealing*|60 °C|30 sec|12|
|Extension|72 °C|30 sec|12|
|Final Extension|72 °C|1 min|1|
|Hold|4 °C | ∞|1|

* Perform a 1X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|Biotin reaction| 25 μl|
|Agencourt® AMPure® XP reagent|  25 μl|
|**Total Volume**| **50 μl**|

* Thoroughly resuspend the beads by pipetting up and down multiple times.
* Incubate the plate/tube at room temperature for 5 – 15 min to allow the DNA to bind to the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Carefully remove and discard 115 μl of supernatant.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol.
* Keeping the plate/tube on the magnet, add 200 μl of 80% ethanol.
* Incubate the plate/tube at room temperature for ≥30 sec.
* Carefully remove and discard the ethanol. Try to remove all residual ethanol without disturbing the beads.
* Dry the beads at room temperature, until all of the ethanol has evaporated. **Caution: over-drying the beads may result in dramatic yield loss.**
* Remove the plate/tube from the magnet.
* Thoroughly resuspend the beads in 25 μl of 10 mM Tris-HCl (pH 8.0). 
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 24 μl of the clear supernatant to a new plate/tube and proceed to next step.


### Quant samples
**Procedure (Standard HS DNA protocol)**
* Set up the required number of 0.5-mL tubes for standards and samples. The Qubit® RNA HS Assay requires 2 standards.
* Label the tube lids.
* Prepare the Qubit® working solution by diluting the Qubit® DNA HS Reagent 1:200 in Qubit® DNA HS Buffer. Use a clean plastic tube each time you prepare Qubit® working solution. **Do not mix the working solution in a glass container.**
* Add 190 μL of Qubit® working solution to each of the tubes used for standards.
* Add 10 μL of each Qubit® standard to the appropriate tube, then mix by vortexing 2–3 seconds. Be careful not to create bubbles.
* Add Qubit® working solution to individual assay tubes so that the final volume in each tube after adding sample is 200 μL.
* Add each sample to the assay tubes containing the correct volume of Qubit® working solution, then mix by vortexing 2–3 seconds. The final volume in each tube should be 200 μL.
* Allow all tubes to incubate at room temperature for 2 minutes.
* On the Home screen of the Qubit® 3.0 Fluorometer, press DNA, then select DNA: High Sensitivity as the assay type. The “Read standards” screen is displayed. Press Read Standards to proceed.
* Insert the tube containing Standard #1 into the sample chamber, close the lid, then press Read standard. When the reading is complete (~3 seconds), remove Standard #1.
* Insert the tube containing Standard #2 into the sample chamber, close the lid, then press Read standard. When the reading is complete, remove Standard #2.
* Press Run samples.
* On the assay screen, select the sample volume and units
* Insert a sample tube into the sample chamber, close the lid, then press Read tube. When the reading is complete (~3 seconds), remove the sample tube.
* Repeat step last step until all samples have been read

### Verify

* Run samples on BioAnalyzer/Tape Station/Fragment analyzer 




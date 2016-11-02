#EecSeq Lab Protocol

The Expressed Exome Capture Sequencing protocol is designed to create exome capture probes directly from RNA.  The probes are then used from hybrid capture of exome DNA sequences, allowing for genotyping of alleles at expressed genes.

[Foo](#foo)

* [RNA Prep](#rna-prep)


###Extract RNA from 2 control and 2 heat shock individuals
####Using unmodified RNAeasy kit (Below are summary steps)
###Quantify all RNA samples
##Visualize RNA of BioAnalyzer
##Begin KAPA Stranded mRNA-Seq Kit ##using 1/2 rxn volumes##
###Additional reagents needed:
###Custom Oligos needed:
### Anneal Adapters
###mRNA Capture
####mRNA Elution, Fragmentation and Priming
###Safe Stopping Point
###1st Strand Synthesis
###2nd Strand Synthesis and Marking
###2nd Strand Synthesis and Marking Cleanup
###SAFE STOPPING POINT
###A-Tailing
###A-Tailing immediately
###A-Tailing after safe stopping point
###Adapter Ligation
####Adapter concentration will vary depending on overall RNA yield, see table below:
####This will be where we insert the custom adapters that are barcoded with RE sites
###Post-Ligation Cleanup
###Safe Stopping Point
###2nd Post-Ligation Cleanup
###SAFE STOPPING POINT
###Library Amplificiation
##NEEDS TO BE TESTED
###Library Amplification Cleanup
## Quant libraries
###Safe Stopping Point
### DSN Normalization
#### DSN needs to be properly dilued and should be tested for activity levels before proceeding
#####This protocol was taken from Illumina's recommendations
###Safe Stopping Point
###SPRI Cleanup
###PCR Enrichment
####Illumina recommends 12 cycles.  This could be increased.
####Alternatively, we could split vials into libraries and probes here
###SPRI Cleanup
## Quant libraries
##Split finished cDNA library for each sample into two vials (8 vials total)
###Safe Stopping Point
##Probe Synthesis
###Remove adapters from cDNA
####Materials needed
###Remove 5' and 3' overhangs
###Safe Stopping Point
##Biotin Labeling
###Materials needed
###Procedure
####Optional:Control reaction
###Safe Stopping Point
##Preparation of whole genome libraries using KAPA HyperPlus Kit
#### This assumes that genomic DNA is already extracted and sheared.
### Anneal Adapters
### End repair
### Adapter ligation
### Post-ligation Cleanup
### Quant samples
###Library Amplification
### Post-amplification Cleanup
###Safe Stopping Point
## Hybridization and Capture
####Materials needed
###Hybridization
### Preparation of Dynabeads
### Washes
### Library re-amplification
### Quant samples
### Verify







##RNA Prep

###Extract RNA from 2 control and 2 heat shock individuals
*Refer to manual during procedure (steps below are for notes and comments)*
####Using unmodified RNAeasy kit (Below are summary steps)
Notes before starting
* If purifying RNA from cell lines rich in RNases, or tissue, add either 10 μl β-mercaptoethanol (β-ME), or 20 μl 2 M dithiothreitol (DTT)*, to 1 ml Buffer RLT. Buffer RLT with β-ME or DTT can be stored at room temperature for up to 1 month.
* Add 4 volumes of ethanol (96–100%) to Buffer RPE for a working solution.
* Remove RNAlater®-stabilized tissue from the reagent using forceps.
Procedure
* Do not use more than 30 mg tissue. Disrupt the tissue and homogenize the lysate in the appropriate volume of Buffer RLT (350 μl for less than 20 mg; 650 μl for less than 30 mg). Centrifuge the lysate for 3 min at maximum speed. Carefully remove the supernatant by pipetting, and use it in step 2.
* Add 1 volume of 70% ethanol to the lysate, and mix well by pipetting. Do not centrifuge. Proceed immediately to step 3.
* Transfer up to 700 μl of the sample, including any precipitate, to an RNeasy Mini spin column placed in a 2 ml collection tube (supplied). Close the lid, and centrifuge for 15 s at ≥8000 x g. Discard the flow-through.
* Add 700 μl Buffer RW1 to the RNeasy spin column. Close the lid, and centrifuge for 15 s at ≥8000 x g. Discard the flow-through.
* Add 500 μl Buffer RPE to the RNeasy spin column. Close the lid, and centrifuge for 15 s at ≥8000 x g. Discard the flow-through.
* Add 500 μl Buffer RPE to the RNeasy spin column. Close the lid, and centrifuge for 2 min at ≥8000 x g.
* **Optional:** Place the RNeasy spin column in a new 2 ml collection tube (supplied). Centrifuge at full speed for 1 min to dry the membrane.
* Place the RNeasy spin column in a new 1.5 ml collection tube (supplied). Add 30–50 μl RNase-free water directly to the spin column membrane. Close the lid, and centrifuge for 1 min at ≥8000 x g to elute the RNA.
* If the expected RNA yield is >30 μg, repeat step 7 using another 30–50 μl of RNase-free water, or using the eluate from step 7 (if high RNA concentration is required). Reuse the collection tube from step 7.

###Quantify all RNA samples
Results will be used for calibration points during library generation
Refer to manual during procedure (steps below are for notes and comments)

**Procedure (Standard HS RNA protocol)**
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

##Visualize RNA of BioAnalyzer

See Becca Certner for latest Vollmer lab protocol

##Begin KAPA Stranded mRNA-Seq Kit ##using 1/2 rxn volumes##
This should take 8-10 hours
Refer to manual during procedure (steps below are for notes and comments)

###Additional reagents needed:

Annealing buffer stock (10X):

| Component| Concentration|
|----------|--------------|
| Tris HCl, pH 8| 100 mM|
| NaCl|500 mM|
| EDTA| 10 mM|
  

###Custom Oligos needed:

|Oligo Name| Sequence|
|----------|---------|
|RNA_P2.1_H3|P*CAAGCTTAGATCGGAAGAGCGAGAACAA
|RNA_P2.1_NC|P*CCCATGGAGATCGGAAGAGCGAGAACAA
|RNA_P2.1_SA|P*CGTCGACAGATCGGAAGAGCGAGAACAA
|RNA_P2.1_BS|P*CTGTACAAGATCGGAAGAGCGAGAACAA
|RNA_P2.2_H3|GTGACTGGAGTTCACACGTGTGCTCTTCCGATCTTTCGAAG*T
|RNA_P2.2_NC|GTGACTGGAGTTCACACGTGTGCTCTTCCGATCTGGTACCG*T
|RNA_P2.2_SA|GTGACTGGAGTTCACACGTGTGCTCTTCCGATCTCAGCTGG*T
|RNA_P2.2_BS|GTGACTGGAGTTCACACGTGTGCTCTTCCGATCTACATGTG*T
|RNA_P1.1_H3|ACACTCTTTCCCTACACGACGCTCTTCCGATCTAAGCTTG*T
|RNA_P1.1_NC|ACACTCTTTCCCTACACGACGCTCTTCCGATCTCCATGGG*T
|RNA_P1.1_SA|ACACTCTTTCCCTACACGACGCTCTTCCGATCTGTCGACG*T
|RNA_P1.1_BS|ACACTCTTTCCCTACACGACGCTCTTCCGATCTTGTACAG*T
|RNA_P1.2_H3|P*CTTCGAAATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
|RNA_P1.2_NC|P*CGGTACCATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
|RNA_P1.2_SA|P*CCAGCTGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
|RNA_P1.2_BS|P*CACATGTATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 

### Anneal Adapters
Single-stranded oligos need to be annealed with their appropriate partner before ligation. 
 
* To create Adapter P1, combine each oligo 1.1 with its complementary oligo 1.2 in a 1:1 ratio in working strength annealing buffer (final buffer concentration 1x) for a total annealed adapter concentration of 40uM (for example, if purchased oligos are resuspended to an initial concentration of 100uM, use 40ul oligo 1.1, 40ul oligo 1.2, 10ul 10x annealing buffer and 10ul nuclease-free water). Do the same for oligos 2.1 and 2.2 to create the common adapter P2. 
2.   In a thermocyler, incubate at 97.5°C for 2.5 minutes, and then cool at a rate of not greater than 3°C per minute until the solution reaches a temperature of 21°C. Hold at 4°C.
3.   Prepare final working strength concentrations of annealed adapters from this annealed stock (the appropriate working stock dilution for your experiment can be determined from our ligation molarity calculator). For convenience, it is possible to store the adapters at 4°C while in active use.


**Procedure:**
* Prepare mastermixes for 4 libraries
  * See tables in manual and tables below for guidelines:

|Component| Total Volume Needed for 4 RXNs|
|---------|--------------------|
|**1st Strand Synthesis Master Mix:**|--|
|1st Strand Synthesis Buffer|22 μl|
|KAPA Script|2 μl|
|**Total Master Mix Volume**| **24 μl**|
|**Final reaction composition:**|--|
|1st Strand Synthesis Master Mix| 5 μl|
|Fragmented, primed RNA|10 μl| 
|**Total Reaction Volume**| **15 μl**|

|Component| Total Volume Needed for 4 RXNs|
|---------|--------------------|
|**2nd Strand Synthesis and Marking Master Mix:**|--|
|2nd Strand Marking Buffer|62 μl|
|2nd Strand Synthesis Enzyme Mix|4 μl|
|**Total Master Mix Volume**| **66 μl**|
|**Final reaction composition:**|--|
|2nd Strand Synthesis and Marking Master Mix| 15 μl|
|Fragmented, primed RNA|15 μl| 
|**Total Reaction Volume**| **30 μl**|

|Component| Total Volume Needed for 4 RXNs (10% excess)|
|---------|--------------------|
|**A-Tailing Master Mix:**|--|
|Water|52.8 μl|
|10X KAPA A-Tailing Buffer|6.6 μl|
|KAPA A-Tailing Enzyme|6.6 μl|
|**Total Master Mix Volume**| **66 μl**|
|**Resuspend beads in a volume of:**|** 15 μl** |

|Component| Total Volume Needed for 4 RXNs|
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

|Component| Total Volume Needed for 4 RXNs (10% excess)|
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

###mRNA Capture
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
####mRNA Elution, Fragmentation and Priming
* Prepare the required volume of 1X Fragment, Prime and Elute Buffer as follows:

|Component|Volume per sample|
|---------|------|
|Water| 5.5 μl|
|Fragment, Prime and Elute Buffer (2X)| 5.5 μl|
|**Total Volume**| **11 μl**|

* Thoroughly resuspend the KAPA mRNA Capture Beads with captured mRNA prepared in Step 2.13 above in 11 μl of 1X Fragment, Prime and Elute Buffer.

---
###Safe Stopping Point
Resuspended beads with captured mRNA may be stored at 4 oC for up to 24 hours. Do not freeze the samples as this will damage the beads. When ready, proceed to step below.

---

* Place the plate/tubes in a thermal cycler and carry out the fragmentation and priming program as follows:

|Desired Fragment Size| Temp.| Duration|
|---------------------|------|---------|
|100 – 200 bp|94 °C|8 min|
|**200 – 300 bp**|**94 °C**|**6 min**|
|300 – 400 bp|85 °C|6 min|

---
**_We are planning on using 125 bp or 150 bp sequencing, so I think it makes sense to use the 200-300 range_**

---

* Immediately place the plate/tube on a magnet to capture the beads, and incubate until the liquid is clear. **Caution: To prevent hybridization of poly(A)- rich RNA to the capture beads, do not allow the sample to cool before placing on the magnet.**
* Carefully remove 10 μl of the supernatant containing the eluted, fragmented, and primed RNA into a separate plate or tube.
* Proceed immediately to **1st Strand Synthesis**.

###1st Strand Synthesis
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

###2nd Strand Synthesis and Marking
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

###2nd Strand Synthesis and Marking Cleanup

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

###SAFE STOPPING POINT
* Resuspend the beads in 15 μl 1X A-Tailing Buffer (Table 5B), cover the reaction and store at 4 oC for up to 24 hours. Do not freeze the samples as this will damage the AMPure® XP® beads. When ready, proceed to **A-Tailing after Safe Stopping Point**.

---

###A-Tailing
* A-Tailing is performed either directly after the 2nd Strand Synthesis and Marking Cleanup, or after the Safe Stopping Point, where beads were resuspended in 1X A-Tailing Buffer and stored at 4 °C for up to 24 hours. Depending on your chosen work ow, proceed with either **A-Tailing immediately** or **Section 7B: A-Tailing after Safe Stopping Point**.

###A-Tailing immediately
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

###A-Tailing after safe stopping point
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

###Adapter Ligation
* Set up the adapter ligation reactions as follows:

|Component|Volume|
|---------|------|
|Beads with A-tailed DNA| 15 μl |
|Adapter Ligation Master Mix| 17.5 μl |
|**Adapters***| 2.5 μl|
|**Total Volume**| **35 μl**|

####Adapter concentration will vary depending on overall RNA yield, see table below:
|Quantity of starting material|Adapter stock concentration|Adapter concentration in ligation reaction|
|----|----|-----|
|100 – 250 ng|140 nM| 10 nM|
|251 – 500 ng|350 nM|25 nM|
|501 – 2000 ng|700 nM|50 nM|
|2001 – 4000 ng|1400 nM|100 nM|

####This will be where we insert the custom adapters that are barcoded with RE sites

* Mix thoroughly by pipetting up and down several times to resuspend the beads.
* Incubate the plate/tube at 20 °C for 15 min.
* Proceed immediately to **1st Post-Ligation Cleanup**.

###Post-Ligation Cleanup

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

###Safe Stopping Point
The solution with resuspended beads can be stored at 4 °C for up to 24 hours. Do not freeze the beads, as this can result in dramatic loss of DNA. When ready, proceed to **2nd Post-Ligation Cleanup**.

---

###2nd Post-Ligation Cleanup

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
###SAFE STOPPING POINT
The purified, adapter-ligated library DNA may be stored at 4 °C for up to 1 week, or frozen at -20 °C for up to 1 month. When ready, proceed to **Library Amplification**.

---
###Library Amplificiation

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
|Denaturation|98 °C|15 sec|16|
|Annealing*|60 °C|30 sec|16|
|Extension|72 °C|30 sec|16|
|Final Extension|72 °C|5 min|1|
|Hold|10 °C | ∞|1|

---
##NEEDS TO BE TESTED
*Annealing temperature may need to be optimized*

I think that 8 cycles should be enough.  It's half the maximum recommended number for the kit, and given the that we have another PCR step after DSN, I don't think we need to maximize here

---

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

###Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---


### DSN Normalization
####DSN needs to be properly dilued and should be tested for activity levels before proceeding

#####This protocol was taken from Illumina's recommendations

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

* Caution
  * Following incubation, keep the thermal cycler lid closed and the temperature held at 68°C. Do not remove the reaction mix tube from thermal cycler prior to and during DSN treatment.
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

###Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---
###SPRI Cleanup

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

###PCR Enrichment

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
---
####Illumina recommends 12 cycles.  This could be increased.
####Alternatively, we could split vials into libraries and probes here
---

###SPRI Cleanup

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
* Thoroughly resuspend the beads in X μl of 10 mM Tris-HCl (pH 8.0).
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 24 μl of the clear supernatant to a new plate/tube and proceed to next step.

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

##Split finished cDNA library for each sample into two vials (8 vials total)
* One tube for sequencing
* One tube for probe synthesis

---

###Safe Stopping Point
This is a safe stopping point. If you are stopping, store your sample at ‐15° to ‐25°C.

---

##Probe Synthesis

Review quantifications for probes.  We would like around 1000-2000 ng of probes per capture.  If not enough of probes is obtained, the PCR product can be re-amplified.

###Remove adapters from cDNA

####Materials needed
| reagent                                | producer         | serial #      |            
|----------------------------------------|------------------|-----------|
|BsrGI-HF| NEB| R3575S|
|Mung Bean Nuclease| NEB| M0250S|
|SalI-HF| NEB| R3138T|
|NcoI-HF| NEB| R3193S|
|HindIII-HF| NEB| R3104T|
|Agencourt AMPure XP  |Beckman Coulter   | A63881|

* Setup four different restriction digests

|Component|Volume||Component|Volume||Component|Volume||Component|Volume|
|---------|------||---------|------||---------|------||---------|------|
| 10X Cutsmart Buffer| 5 μl| | 10X Cutsmart Buffer| 5 μl| | 10X Cutsmart Buffer| 5 μl| | 10X Cutsmart Buffer| 5 μl| 
| BsrGI-HF Enzyme (100 units)| 5 μl || SalI-HF Enzyme (100 units)| 1 μl | | NcoI-HF Enzyme (100 units)| 5 μl | | HindIII-HF Enzyme (100 units)| 1 μl|  
| Molecular Grade H20| 16 μl|| Molecular Grade H20| 21 μl|| Molecular Grade H20| 16 μl|| Molecular Grade H20| 21 μl|
| DSN Enriched Library| 24 μl| | DSN Enriched Library| 24 μl| | DSN Enriched Library| 24 μl| | DSN Enriched Library| 24 μl| 

**Make sure library is matched with appopriate enzyme**

* Incubate reactions in thermocycler at 37°C for at least 8 hours

* Perform a 1.5X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|RE treated DSN Library| 50 μl|
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
* Thoroughly resuspend the beads in 20 μl of 10 mM Tris-HCl (pH 8.0). 
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.

###Remove 5' and 3' overhangs

* **Using the same tubes from the previous step** combine the following:

|Component|Volume|
|---------|------|
|Library| 20 μl|
|10X Mung Bean Nuclease buffer| 5 μl|
|Mung Bean Nuclease (10 units per μl)| 5 μl|
|Molecular Grade H20| 20 μl|
|**Total Volume**| **50 μl**|

* Incubate at 30°C for 30 minutes

* Perform a 1.5X SPRI® cleanup by combining the following:

|Component|Volume|
|---------|------|
|MBN reaction| 50 μl|
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
* Thoroughly resuspend the beads in 21 μl of 10 mM Tris-HCl (pH 8.0). Volume needed depends on the number of captures. Calculate 10 μl per capture plus an aliquot for checking the probes concentration using Qubit.
* Incubate the plate/tube at room temperature for 2 min to allow the DNA to elute off the beads.
* Place the plate/tube on a magnet to capture the beads. Incubate until the liquid is clear.
* Transfer 20 μl of the clear supernatant to a new plate/tube and proceed to next step.

---

###Safe Stopping Point
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
* Replace with 30 μl of H20 
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




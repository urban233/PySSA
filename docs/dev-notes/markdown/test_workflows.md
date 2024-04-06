# Workflow Results
* 
  * [Test workflow 1](#Test-workflow-1) Result: OK (05.04.2024)
  * [Test workflow 2](#Test-workflow-2) Result: OK (04.04.2024; no more testing needed)
  * [Test workflow 3](#Test-workflow-3) Result: OK (05.04.2024)
  * [Test workflow 4](#Test-workflow-4) Result: OK (05.04.2024)
  * [Test workflow 5](#Test-workflow-5) Result: OK (05.04.2024)
* Proteins
  * Import
    * [Test workflow 6](#Test-workflow-6) Result: OK (05.04.2024)
    * [Test workflow 7](#Test-workflow-7) Result: OK (05.04.2024)
    * [Test workflow 8](#Test-workflow-8) Result: OK (05.04.2024)
    * [Test workflow 9](#Test-workflow-9) Result: OK (05.04.2024)
* Sequences
  * Add
    * [Test workflow 10](#Test-workflow-10) Result: OK (06.04.2024)
  * Import
    * [Test workflow 11](#Test-workflow-11) Result: OK (06.04.2024)
* Prediction
  * Monomer
    * [Test workflow 12](#Test-workflow-12) Result: **Error** -> Project menu option is disabled
    * [Test workflow 13](#Test-workflow-13) Result: 
    * [Test workflow 14](#Test-workflow-14) Result: **Error** -> Project menu option is disabled
    * [Test workflow 15](#Test-workflow-15) Result:
  * Multimer
    * [Test workflow 16](#Test-workflow-16) Result:
    * [Test workflow 17](#Test-workflow-17) Result:
    * [Test workflow 18](#Test-workflow-18) Result:
    * [Test workflow 19](#Test-workflow-19) Result:
* Prediction + Analysis
  * Monomer
    * [Test workflow 20](#Test-workflow-20) Result: **Error** -> Project menu option is disabled
    * [Test workflow 21](#Test-workflow-21) Result:
  * Multimer
    * [Test workflow 22](#Test-workflow-22) Result:
    * [Test workflow 23](#Test-workflow-23) Result:
* Distance Analysis
  * [Test workflow 24](#Test-workflow-24) Result:
  * [Test workflow 25](#Test-workflow-25) Result:

# Test Workflows (Extended)
### Test workflow 1
 * Start PySSA
 * Help/Documentation
   * Error: No focus, no windows activation

### Test workflow 3
 * Start PySSA
 * Empty Workspace
 * Only New and Import
   * Error: Delete etc. is visible

### Test workflow 4
 * Start PySSA
 * New project: "Kalata-B1-Prediction"
 * Add protein sequence: CGETCVGGTCNTPGCTCSWPVCTRNGLPV
   * Error: All menue entries become visible.

### Test workflow 5
 * Start PySSA
 * New project: "Kalata-B1-Prediction"
 * Add protein sequence: CGETCVGGTCNTPGCTCSWPVCTRNGLPV
 * Menu: Prediction becomes activated
 * Prediction/Monomer: Selection of sequence(s)
   * Error: Function is disconnected from sequence entries

# Proteins
## Import
### Test workflow 6
Description: This workflow outlines the steps for adding a protein with a PDB id
to an existing project.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * A test project with the name "F_BMP2-xray_mono" exists.
 * The test project does not contain a protein with the name "3BMP".

Steps:
 1. Open project "F_BMP2-xray_mono".
 2. Switch tab to "Proteins".
 3. Import protein with the id "3bmp".

Post-conditions:
 * Project contains a protein with the name 3BMP (case-sensitive!).

Cleanup:
 * Delete the protein "3BMP".
 * Exit application.

Post-conditions after cleanup:
 * Project has no proteins (tree view on proteins tab is empty).


### Test workflow 7

Description: This workflow outlines the steps for adding a protein with a PDB id
to a new project.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name "BMP2-xray_mono".

Steps:
 1. Create project "BMP2-xray_mono".
 2. Switch tab to "Proteins".
 3. Import protein with the id "3bmp".

Post-conditions:
 * Project contains a protein with the name 3BMP (case-sensitive!).

Cleanup:
 * Close the project.
 * Delete the project "BMP2-xray_mono".
 * Exit application.

Post-conditions after cleanup:
 * Project "BMP2-xray_mono" is removed from the workspace.


### Test workflow 8

Description: This workflow outlines the steps for adding a protein with an
existing .pdb file to an existing project.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * A test project with the name "F_BMP2-xray_mono" exists.
 * The test project does not contain a protein with the name "3BMP".

Steps:
 1. Open project "BMP2-xray_mono".
 2. Switch tab to "Proteins".
 3. Import protein from the "3bmp.pdb" file.

Post-conditions:
 * Project contains a protein with the name 3bmp (case-sensitive!).

Cleanup:
 * Delete the protein "3bmp".
 * Exit application.

Post-conditions after cleanup:
 * Project has no proteins (tree view on proteins tab is empty).


### Test workflow 9

Description: This workflow outlines the steps for adding a protein with an
existing .pdb file to a new project.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name "BMP2-xray_mono".

Steps:
 1. Create project "BMP2-xray_mono".
 2. Switch tab to "Proteins".
 3. Import protein from the "3bmp.pdb" file.

Post-conditions:
 * Project contains a protein with the name 3bmp (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "BMP2-xray_mono".
 * Exit application.

Post-conditions after cleanup:
 * Project "BMP2-xray_mono" is removed from the workspace.

# Sequences
## Add
### Test workflow 10
Description: This workflow outlines the steps for adding a sequence with an
existing sequence string to a new project.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name "Kalata-B1-Seqs".

Steps:
 1. Create project "Kalata-B1-Seqs".
 2. Click on the "add sequence icon".
 3. Enter the sequence name "kalata_b1" and click on "Next".
 4. Paste the sequence "CGETCVGGTCNTPGCTCSWPVCTRNGLPV " into the text field.
 5. Click on "Add"

Post-conditions:
 * Project contains a sequence with the name kalata_b1 (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Seqs".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Seqs" is removed from the workspace.


## Import
### Test workflow 11
Description: This workflow outlines the steps for adding a sequence with an
existing .fasta file to a new project.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name "Kalata-B1-Seqs".

Steps:
 1. Create project "Kalata-B1-Seqs".
 2. Click on the "import sequence icon".
 3. Choose the kalata.fasta file from the filesystem.
 4. Click on "Import"

Post-conditions:
 * Project contains a sequence with the name Kalata (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Seqs".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Seqs" is removed from the workspace.

# Prediction
## Monomer
### Test workflow 12
Description: This workflow outlines the steps for predicting a monomeric protein
structure based on a new project where the sequence is added through the "Add"
functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-B1-Prediction".

Steps:
 1. Create project "Kalata-B1-Prediction".
 2. Add sequence "kalata_b1" with the sequence "CGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 3. Click on the newly added sequence "kalata_b1".
 4. Click on "Monomer" in the menu.
 5. Click on "Predict".
 6. Wait for the prediction to finish!

Post-conditions:
 * Project contains a protein with the name kalata_b1 (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Prediction".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Prediction" is removed from the workspace.


### Test workflow 13
Description: This workflow outlines the steps for predicting a monomeric protein
structure based on a new project where the sequence is added through the
"Import" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-B1-Prediction".

Steps:
 1. Create project "Kalata-B1-Prediction".
 2. Import sequence "kalata" with the kalata.fasta file.
 3. Click on the newly added sequence "kalata".
 4. Click on "Monomer" in the menu.
 5. Click on "Predict".
 6. Wait for the prediction to finish!

Post-conditions:
 * Project contains a protein with the name kalata (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Prediction".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Prediction" is removed from the workspace.


### Test workflow 14
Description: This workflow outlines the steps for predicting two monomeric
protein structures based on a new project where the sequences are added through
the "Add" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-B1-Prediction".

Steps:
 1. Create project "Kalata-B1-Prediction".
 2. Add sequence "kalata_b1" with the sequence "CGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 3. Add sequence "kalata_b1_m" with the sequence
    "GGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 4. Click on "Monomer" in the menu.
 5. Click on "Predict".
 6. Wait for the prediction to finish!

Post-conditions:
 * Project contains two proteins with the name kalata_b1 and kalata_b1_m
   (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Prediction".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Prediction" is removed from the workspace.


### Test workflow 15
Description: This workflow outlines the steps for predicting two monomeric
protein structures based on a new project where the sequences are added through
the "Import" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-B1-Prediction".

Steps:
 1. Create project "Kalata-B1-Prediction".
 2. Import sequence "kalata" with the kalata.fasta file.
 3. Import sequence "kalata_m" with the kalata_m.fasta file.
 4. Click on "Monomer" in the menu.
 5. Click on "Predict".
 6. Wait for the prediction to finish!

Post-conditions:
 * Project contains a protein with the name kalata (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Prediction".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Prediction" is removed from the workspace.


## Multimer
### Test workflow 16
Description: This workflow outlines the steps for predicting a dimeric protein
structure based on a new project where the sequence is added through the "Add"
functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-Dimer-Prediction".

Steps:
 1. Create project "Kalata-Dimer-Prediction".
 2. Add sequence "kalata_b1_dimer" with the sequence
    "CGETCVGGTCNTPGCTCSWPVCTRNGLPV,CGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 3. Click on the newly added sequence "kalata_b1_dimer".
 4. Click on "Multimer" in the menu.
 5. Click on "Predict".
 6. Wait for the prediction to finish!

Post-conditions:
 * Project contains a protein with the name kalata_b1_dimer (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-Dimer-Prediction".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-Dimer-Prediction" is removed from the workspace.


### Test workflow 17
Description: This workflow outlines the steps for predicting a dimeric protein
structure based on a new project where the sequence is added through the
"Import" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-Dimer-Prediction".

Steps:
 1. Create project "Kalata-Dimer-Prediction".
 2. Import sequence "kalata_dimer" with the kalata_dimer.fasta file.
 3. Click on the newly added sequence "kalata_dimer".
 4. Click on "Multimer" in the menu.
 5. Click on "Predict".
 6. Wait for the prediction to finish!

Post-conditions:
 * Project contains a protein with the name kalata_dimer (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-Dimer-Prediction".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-Dimer-Prediction" is removed from the workspace.


### Test workflow 18
Description: This workflow outlines the steps for predicting two dimeric protein
structures based on a new project where the sequences are added through the
"Add" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-Dimer-Prediction".

Steps:
 1. Create project "Kalata-Dimer-Prediction".
 2. Add sequence "kalata_b1_dimer" with the sequence
    "CGETCVGGTCNTPGCTCSWPVCTRNGLPV,CGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 3. Add sequence "kalata_b1_dimer_m" with the sequence
    "GGETCVGGTCNTPGCTCSWPVCTRNGLPV,GGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 4. Click on "Multimer" in the menu.
 5. Click on "Predict".
 6. Wait for the prediction to finish!

Post-conditions:
 * Project contains two proteins with the name kalata_b1 and kalata_b1_m
   (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-Dimer-Prediction".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-Dimer-Prediction" is removed from the workspace.


### Test workflow 19
Description: This workflow outlines the steps for predicting two dimeric protein
structures based on a new project where the sequences are added through the
"Import" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-Dimer-Prediction".

Steps:
 1. Create project "Kalata-Dimer-Prediction".
 2. Import sequence "kalata_dimer" with the kalata_dimer.fasta file.
 3. Import sequence "kalata_dimer_m" with the kalata_dimer_m.fasta file.
 4. Click on "Multimer" in the menu.
 5. Click on "Predict".
 6. Wait for the prediction to finish!

Post-conditions:
 * Project contains a protein with the name kalata (case-sensitive!).



Cleanup:
 * Close the project.
 * Delete the project "Kalata-Dimer-Prediction".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-Dimer-Prediction" is removed from the workspace.

# Prediction + Analysis
## Monomer
### Test workflow 20
Description: This workflow outlines the steps for predicting and analyzing a
monomeric protein structure based on a new project where the sequence is added
through the "Add" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-B1-Prediction_Analysis".

Steps:
 1. Create project "Kalata-B1-Prediction_Analysis".
 2. Add sequence "kalata_b1" with the sequence "CGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 3. Switch tab to "Proteins".
 4. Import protein from the "1NB1" PDB id.
 5. Click on "Monomer" in the menu.
 6. Tick checkbox "Add Analysis".
 7. Setup analysis run (kalata_b1 and 1NB1)
 8. Start prediction + analysis job.
 9. Wait for the job to finish!

Post-conditions:
 * Project contains a protein with the name kalata_b1 (case-sensitive!).
 * Project contains a protein pair with the name kalata_b1_A_vs_1NB1_A



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Prediction_Analysis".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Prediction_Analysis" is removed from the workspace.


### Test workflow 21
Description: This workflow outlines the steps for predicting and analyzing two
monomeric protein structures based on a new project where the sequence is added
through the "Add" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-B1-Prediction_Analysis".

Steps:
 1.  Create project "Kalata-B1-Prediction_Analysis".
 2.  Add sequence "kalata_b1" with the sequence "CGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 3.  Add sequence "kalata_b1_m" with the sequence
     "GGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 4.  Switch tab to "Proteins".
 5.  Import protein from the "1NB1" PDB id.
 6.  Click on "Monomer" in the menu.
 7.  Tick checkbox "Add Analysis".
 8.  Setup analysis run (kalata_b1 and 1NB1)
 9.  Start prediction + analysis job.
 10. Wait for the job to finish!

Post-conditions:
 * Project contains two new proteins with the name kalata_b1 and kalata_b1_m
   (case-sensitive!)
 * Project contains a protein pair with the name kalata_b1_A_vs_1NB1_A



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Prediction_Analysis".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Prediction_Analysis" is removed from the workspace.


## Multimer
### Test workflow 22
Description: This workflow outlines the steps for predicting and analyzing a
dimeric protein structure based on a new project where the sequence is added
through the "Add" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-B1-Prediction_Analysis".

Steps:
 1. Create project "Kalata-B1-Prediction_Analysis".
 2. Add sequence "kalata_b1_dimer" with the sequence
    "CGETCVGGTCNTPGCTCSWPVCTRNGLPV,CGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 3. Switch tab to "Proteins".
 4. Import protein from the "1NB1" PDB id.
 5. Click on "Multimer" in the menu.
 6. Tick checkbox "Add Analysis".
 7. Setup analysis run (kalata_b1_dimer and 1NB1)
 8. Start prediction + analysis job.
 9. Wait for the job to finish!

Post-conditions:
 * Project contains a protein with the name kalata_b1 (case-sensitive!).
 * Project contains a protein pair with the name kalata_b1_dimer_A_vs_1NB1_A



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Prediction_Analysis".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Prediction_Analysis" is removed from the workspace.


### Test workflow 23
Description: This workflow outlines the steps for predicting and analyzing two
dimeric protein structures based on a new project where the sequence is added
through the "Add" functionality.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name
   "Kalata-B1-Prediction_Analysis".

Steps:
 1.  Create project "Kalata-B1-Prediction_Analysis".
 2.  Add sequence "kalata_b1_dimer" with the sequence
     "CGETCVGGTCNTPGCTCSWPVCTRNGLPV,CGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 3.  Add sequence "kalata_b1_dimer_m" with the sequence
     "GGETCVGGTCNTPGCTCSWPVCTRNGLPV,GGETCVGGTCNTPGCTCSWPVCTRNGLPV".
 4.  Switch tab to "Proteins".
 5.  Import protein from the "1NB1" PDB id.
 6.  Click on "Multimer" in the menu.
 7.  Tick checkbox "Add Analysis".
 8.  Setup analysis run (kalata_b1_dimer and 1NB1)
 9.  Start prediction + analysis job.
 10. Wait for the job to finish!

Post-conditions:
 * Project contains two new proteins with the name kalata_b1 and kalata_b1_m
   (case-sensitive!)
 * Project contains a protein pair with the name kalata_b1_dimer_A_vs_1NB1_A



Cleanup:
 * Close the project.
 * Delete the project "Kalata-B1-Prediction_Analysis".
 * Exit application.

Post-conditions after cleanup:
 * Project "Kalata-B1-Prediction_Analysis" is removed from the workspace.

# Distance Analysis
### Test workflow 24
Description: This workflow outlines the steps for analyzing two protein
structures based on a new project where two proteins are imported with their PDB
id.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name "BMP2-xray_analysis".

Steps:
 1. Create project "BMP2-xray_analysis".
 2. Switch tab to "Proteins".
 3. Import protein with the id "3bmp".
 4. Import protein with the id "6omn".
 5. Click on "Distance" in the menu.
 6. Setup analysis run (3BMP and 6OMN)
 7. Start analysis job.
 8. Wait for the job to finish!

Post-conditions:
 * Project contains a protein pair with the name 3BMP_A_vs_6OMN_E



Cleanup:
 * Close the project.
 * Delete the project "BMP2-xray_analysis".
 * Exit application.

Post-conditions after cleanup:
 * Project "BMP2-xray_analysis" is removed from the workspace.


### Test workflow 25
Description: This workflow outlines the steps for analyzing two protein
structures (creating two analysis runs) based on a new project where two
proteins are imported with their PDB id.

Preconditions:
 * PySSA and PyMOL are started and windows are arranged.
 * The workspace does not contain a project with the name "BMP2-xray_analysis".

Steps:
 1. Create project "BMP2-xray_analysis".
 2. Switch tab to "Proteins".
 3. Import protein with the id "3bmp".
 4. Import protein with the id "6omn".
 5. Click on "Distance" in the menu.
 6. Setup analysis run (3BMP and 6OMN)
 7. Setup analysis run (3BMP (A) and 6OMN (F))
 8. Start analysis job.
 9. Wait for the job to finish!

Post-conditions:
 * Project contains a protein pair with the name 3BMP_A_vs_6OMN_F



Cleanup:
 * Close the project.
 * Delete the project "BMP2-xray_analysis".
 * Exit application.

Post-conditions after cleanup:
 * Project "BMP2-xray_analysis" is removed from the workspace.



--------------------------------------------------------------------------------
### Test workflow 2
 * Start PySSA
 * Change windows
 * Rearrange windows
OK



----------
Tags: 
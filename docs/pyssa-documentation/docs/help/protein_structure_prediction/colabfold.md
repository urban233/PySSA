??? general-information "General Information"
    
    Here are general information about the _ColabFold Prediction_ dialog window.

# ColabFold Prediction
<div class="grid cards" markdown>

-   __Overview__

     In this dialog you can ...

    - predict monomeric or/and multimeric protein structures.
    - predict protein structures with subsequent distance analysis'.

</div>

---
## Getting Started
### Run a protein structure prediction
1. [Add](../sequences/sequence_add.md) or [import](../sequences/sequence_import.md) at least one sequence.
2. Under the _Prediction_ menu, click on **Monomer** or **Multimer** depending on the type of sequence.
3. Click on **Predict**.

### Run a protein structure prediction with a subsequent distance analysis
1. [Add](../sequences/sequence_add.md) or [import](../sequences/sequence_import.md) at least one sequence.
2. [Import](../proteins/protein_import.md) a protein structure that you want to compare to the predicted structure.
3. Under the _Prediction_ menu, click on **Monomer** or **Multimer** depending on the type of sequence.
4. Tick the checkbox besides _Add Analysis_.
5. Click on **Go**, to set up the analysis run.
6. Click on **Add**.
7. Choose a first protein structure from the tree.
8. Click on **Next**.
9. Select another protein structure from the tree.
10. Click on **Add**.
11. Click on **Start**.


## Details
- During the structure prediction it is possible to work as normal.
- It is not possible to run another prediction or analysis while a structure prediction is running.
- The prediction runtime depends on
      - sequence length
      - number of chains
      - number of proteins to predict
      - how much RAM can be used
      - how fast the CPU is
- The prediction cannot use an external graphics card. It only runs on CPU.
- The prediction runtime can vary greatly from job to job depending on the job's setup.


##### See Also
[Add Sequence](../sequences/sequence_add.md) :octicons-square-fill-16: [Import Sequence](../sequences/sequence_import.md) :octicons-square-fill-16: [Import Protein](protein_import.md) :octicons-square-fill-16: [Save Protein](protein_save.md)

---

##### Related Overview
:octicons-square-fill-16: [Sequences](../sequences/index.md) :octicons-square-fill-16: [Proteins](../proteins/index.md) :octicons-square-fill-16: [Protein Pairs](../protein_pairs/index.md)
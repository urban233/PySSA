??? general-information "General Information"
    
    Here are general information about the _Advanced Prediction Configuration_ dialog window.

# Advanced Prediction Configuration
<div class="grid cards" markdown>

-   __Overview__

     In this dialog you can ...

    - turn on/off the AMBER force field.
    - choose a pdb template.

</div>

---
## Details
### AMBER Force Field
- ColabFold can use an AMBER force field for energy minimization.
- The default is that the prediction runs **with** AMBER force field.
- If the prediction fails due to an unknown error, it may be helpful to turn off the force field.

### PDB Template
- ColabFold can use pdb files as template.
- The default is that ColabFold uses the pdb70 database.
- If the prediction fails due to an unknown error, it may be helpful to switch _None_.

##### See Also
[Add Sequence](../sequences/sequence_add.md) :octicons-square-fill-16: [Import Sequence](../sequences/sequence_import.md) :octicons-square-fill-16: [Import Protein](protein_import.md) :octicons-square-fill-16: [Save Protein](protein_save.md)

---

##### Related Overview
:octicons-square-fill-16: [Sequences](../sequences/index.md) :octicons-square-fill-16: [Proteins](../proteins/index.md) :octicons-square-fill-16: [Protein Pairs](../protein_pairs/index.md)
??? general-information "General Information"
    
    Here are general information about the _Settings_ dialog window.

# PySSA Settings
<div class="grid cards" markdown>

-   __Overview__

     In this dialog you can ...

    - change the workspace.
    - enable/disable the automatic startup of the _Help Center_ window.
    - change the _Cycles_ value for the distance analysis.
    - change the _Cutoff_ value for the distance analysis.
    - change the background color.
    - change the ray-tracing renderer.
    - change the ray-trace-mode.
    - change the ray texture.

</div>

---
## Details
### Workspace
- All projects are stored in this folder.
- The default value is C:\Users\<username>\.pyssa\default_workspace.

### Automatic startup of the Help Center window
- Tick the checkbox, if the _Help Center_ should open at application startup.
- The checkbox is ticked by default.

### Cycles value (Distance Analysis)
- Defines the number of outlier rejection cycles for the structure alignment.
- The default value is _0_ cycles, that means all residues are used.

### Cutoff value (Distance Analysis)
- Every distance that is above the cutoff will be removed from the result if the number of cycles is greater 0.
- The default value is _1.0_ but is not used in the distance analysis because the number of cycles are set to 0.

### Background color
- Changes the background color of the current PyMOL session.
- _Black_ is the default value.
- Two colors (black and white) can be used.

### Ray-tracing renderer
- Defines the renderer that is used during the ray-tracing process
- The default renderer is _PyMOL internal renderer_.

### Ray-trace-mode.
- Defines how the colors look in the ray-traced image.
- The default value is _normal color + black outline_
- Experiment with different modes and see what fits best for your use case.

### Ray texture.
- Defines the texture PyMOL uses for the proteins during ray-tracing.
- The default value is _None_.
- Experiment with different textures and see what fits best for your use case.

##### See Also
[Advanded Prediction Configuration](../protein_structure_prediction/advanced_prediction_configurations.md)

---

##### Related Overview
:octicons-square-fill-16: [Settings](index.md)
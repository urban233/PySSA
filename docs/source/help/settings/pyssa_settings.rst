PySSA Settings
==============

.. card:: Description

     In this dialog you can ...

    - change the workspace.
    - change the *Cycles* value for the distance analysis.
    - change the *Cutoff* value for the distance analysis.
    - change the background color.
    - change the ray-tracing renderer.
    - change the ray-trace-mode.
    - change the ray texture.

Details
-------
Workspace
************************************************
- All projects are stored in this folder.
- The default value is C:\Users\<username>\.pyssa\default_workspace.

Cycles value (Distance Analysis)
************************************************
- Defines the number of outlier rejection cycles for the structure alignment.
- The default value is *0* cycles, that means all residues are used.

Cutoff value (Distance Analysis)
************************************************
- Every distance that is above the cutoff will be removed from the result if the number of cycles is greater 0.
- The default value is *1.0* but is not used in the distance analysis because the number of cycles are set to 0.

Background color
************************************************
- Changes the background color of the current PyMOL session.
- *Black* is the default value.
- Two colors (black and white) can be used.

Ray-tracing renderer
************************************************
- Defines the renderer that is used during the ray-tracing process
- The default renderer is *PyMOL internal renderer*.

Ray-trace-mode
************************************************
- Defines how the colors look in the ray-traced image.
- The default value is *normal color + black outline*
- Experiment with different modes and see what fits best for your use case.

Ray texture
************************************************
- Defines the texture PyMOL uses for the proteins during ray-tracing.
- The default value is *None*.
- Experiment with different textures and see what fits best for your use case.

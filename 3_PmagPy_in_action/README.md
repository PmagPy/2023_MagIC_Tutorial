## PmagPy in action on MagIC data

In the PmagPy in action portion of the tutorial we will first work through two examples.

### IODP magnetostratigraphy

[IODP Site U1361 example](magnetostratigraphy/magnetostratigraphy.ipynb) The first example is focused on magnetostratigraphic data from this IODP study: 

> *Tauxe, L., Sugisaki, S., Jiménez-Espejo, F., Escutia, C., Cook, C. P., van de Flierdt, T., & Iwai, M. (2015). Geology of the Wilkes Land sub-basin and stability of the East Antarctic Ice Sheet: Insights from rock magnetism at IODP Site U1361. Earth and Planetary Science Letters, 412, 61– 69. https://doi.org/10.1016/j.epsl.2014.12.034*
    
In this example, we will accomplish the following:
- import data from MagIC format
- extract measurements at a specific demagnetization level
- plot specimen interpretations
- visualize magnetostratigraphy with the GPTS
- analyze and plot anisotropy

### Developing an inclination-flattening corrected paleomagnetic pole

[Nonesuch Formation example](magnetostratigraphy/magnetostratigraphy.ipynb) The second example is focused on paleomagnetic data from fine-grained siliciclastic sedimentary rocks of the Mesoproterozoic (ca. 1075 Ma) Nonesuch Formation associated with the forthcoming paper:

> Slotznick, S.P., Swanson-Hysell, N.L., Zhang, Y.,  Clayton, K., Wellman, C.H., Tosca, N.J., and Strother, P.K. (in revision), Reconstructing the paleoenvironment of an oxygenated Mesoproterozoic shoreline and its record of life, GSA Bulletin.
    
In this example, we will accomplish the following:
- import data from MagIC format
- unpack and filter by both component and tilt ammount
- plot directions
- calculate poles
- conduct common mean tests
- conduct fold tests
- conduct E/I inclination shallowing analysis
- implement "classic" inclination shallowing pole calculation
- implement inclination shallowing uncertainty into paleomagnetic pole

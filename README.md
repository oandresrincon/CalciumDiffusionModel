# CalciumDiffusionModel

# Ca_Dynamics Programs

This repository contains a collection of programs related to calcium dynamics modeling in different fiber types. The following sections provide information on how to run the programs and details about each program's functionality and outputs.

## Installation

Before running the programs, make sure to install the required Python libraries by running the following command in your shell:

```bash
pip install -r requirements.txt
```

List of Programs
Program 1: cts.py

Description: Set constants and calculate the concentration of binding sites at rest.

Outputs:

outputs\cts.xls: Constants used in simulations of fiber types I, IIA, IIX, and IIB.
Program 2: f2c.py

Description: Load fluorescence signals and convert them to [Ca+2].

Inputs:

inputs\f1.xls: Fluorescence in single AP fiber type I.
inputs\f2a.xls: Fluorescence in single AP fiber type IIA.
inputs\f2x.xls: Fluorescence in single AP fiber type IIX.
inputs\f2b.xls: Fluorescence in single AP fiber type IIB.
inputs\f1t.xls: Fluorescence in tetanic AP fiber type I.
inputs\f2bt.xls: Fluorescence in tetanic AP fiber type IIB.
Outputs:

figures\cae.pdf: Figure of experimental [Ca+2] during a single AP.
outputs\cae.xls: Data of experimental [Ca+2] during a single AP.
figures\caet.pdf: Figure of experimental [Ca+2] during tetanic APs in a fiber.
outputs\caet.xls: Data of experimental [Ca+2] during tetanic APs in a fiber.
Program 3: eqs.py

Description: Equations used in the model.

Equations:

Eq_Release (Release of Ca+2 equations)
Eqs_ATP (Adenosine triphosphate equations)
Eqs_Tn (Troponin equations)
Eqs_Dye (Dye equations)
Eqs_Pv (Parvalbumin equations)
Eqs_SERCA (Sarco/endoplasmic reticulum Ca2+-ATPase equations)
Eqs_B (Mitochondrial buffers equations)
Eqs_CSQ (Calsequestrin equations)
Eqs_NCX (Na+/Ca2+ exchanger equations)
Eqs_SOCE (Store-operated Ca2+ entry equations)
Eqs_NCE (MITO Na+/Ca2+ exchanger equations)
Eqs_MCU (MITO Ca2+ uniporter equations)
Eqs_MITO (MITO equations, including Eqs_B, Eqs_NCE, and Eqs_MCU)
Eqs_Diff (Diffusion equations)
Program 4: SCMI.py

Description: Generate the figures and data obtained with the single compartment model in fiber type I.

Outputs:

outputs\caatps1.pdf: Figure of [CaATP] during a single AP in fiber type I with the SCM.
outputs\catns1.pdf: Figure of [CaTn] during a single AP in fiber type I with the SCM.
outputs\cadyes1.pdf: Figure of [CaDye] during a single AP in fiber type I with the SCM.
outputs\capvs1.pdf: Figure of [CaPv] during a single AP in fiber type I with the SCM.
outputs\casercas1.pdf: Figure of [CaSERCA] during a single AP in fiber type I with the SCM.
outputs\cancxs1.pdf: Figure of [CaNCX] during a single AP in fiber type I with the SCM.
outputs\camitos1.pdf: Figure of [Ca2+]mito during a single AP in fiber type I with the SCM.
outputs\cats1.pdf: Figure of [Ca2+]total during a single AP in fiber type I with the SCM.
outputs\jrel1.pdf: Figure of release rate of Ca2+ during a single AP in fiber type I with the SCM.
outputs\cs1.xls: Table with the data that generate previous figures.
Program 5: SCMIIA.py

Description: Generate the figures and data obtained with the single compartment model in fiber type IIA.

Outputs:

outputs\caatps2a.pdf: Figure of [CaATP] during a single AP in fiber type IIA with the SCM.
outputs\catns2a.pdf: Figure of [CaTn] during a single AP in fiber type IIA with the SCM.
outputs\cadyes2a.pdf: Figure of [CaDye] during a single AP in fiber type IIA with the SCM.
outputs\capvs2a.pdf: Figure of [CaPv] during a single AP in fiber type IIA with the SCM.
outputs\casercas2a.pdf: Figure of [CaSERCA] during a single AP in fiber type IIA with the SCM.
outputs\cancxs2a.pdf: Figure of [CaNCX] during a single AP in fiber type IIA with the SCM.
outputs\camitos2a.pdf: Figure of [Ca2+]mito during a single AP in fiber type IIA with the SCM.
outputs\cats2a.pdf: Figure of [Ca2+]total during a single AP in fiber type IIA with the SCM.
outputs\jrel2a.pdf: Figure of release rate of Ca2+ during a single AP in fiber type IIA with the SCM.
outputs\cs2a.xls: Table with the data that generate previous figures.
Program 6: SCMIIX.py

Description: Generate the figures and data obtained with the single compartment model in fiber type IIX.

Outputs:

outputs\caatps2x.pdf: Figure of [CaATP] during a single AP in fiber type IIX with the SCM.
outputs\catns2x.pdf: Figure of [CaTn] during a single AP in fiber type IIX with the SCM.
outputs\cadyes2x.pdf: Figure of [CaDye] during a single AP in fiber type IIX with the SCM.
outputs\capvs2x.pdf: Figure of [CaPv] during a single AP in fiber type IIX with the SCM.
outputs\casercas2x.pdf: Figure of [CaSERCA] during a single AP in fiber type IIX with the SCM.
`outputs\cancxs2x.pdf


# CalciumDiffusionModel

#-----------------------------------------------------------------------------
Information to run programs in folder Ca_Dynamics
#-----------------------------------------------------------------------------
-File requierements.txt contain the python libraries used in the folder.
To install: cd to the directory where requirements.txt is located and run: 
pip install -r requirements.txt in shell. 
#-----------------------------------------------------------------------------
-List of programs contained in the folder:
#-----------------------------------------------------------------------------
 Program 1. cts.py

-Description:set constants and calculate concentration of binding sites at
 rest.

-Outputs:
 -outputs\cts.xls (constants used in simulations of fiber types I, IIA, IIX
  and IIB)
#-----------------------------------------------------------------------------
 Program 2. f2c.py

-Description:load fluorescence signals and converts them to [Ca+2]

-Inputs:
 -inputs\f1.xls (fluorescence in single AP fiber type I)
 -inputs\f2a.xls (fluorescence in single AP fiber type IIA)
 -inputs\f2x.xls (fluorescence in single AP fiber type IIX)
 -inputs\f2b.xls (fluorescence in single AP fiber type IIB)
 -inputs\f1t.xls (fluorescence in tetanic AP fiber tipe I)
 -inputs\f2bt.xls (fluorescence in tetanic AP fiber tipe IIB)		
-Outputs:
 -figures\cae.pdf (figure of experimental [Ca+2] during single AP)
 -outputs\cae.xls (data of experimental [Ca+2] during single AP)
 -figures\caet.pdf (figure of experimental [Ca+2] during tetanic APs in fiber)
 -outputs\caet.xls (data of experimental [Ca+2] during tetanic APs in fiber)
#-----------------------------------------------------------------------------
 Program 3. eqs.py
#-----------------------------------------------------------------------------
-Description:equations used in the model.

-Eq_Release (Release of Ca+2 equations)
-Eqs_ATP: (adenosine triphosphate equations)
-Eqs_Tn: (Troponin equations)
-Eqs_Dye: (Dye equations)
-Eqs_Pv: (Parvalbumin equations)
-Eqs_SERCA: (sarco/endoplasmic reticulum Ca2+-ATPase equations)
-Eqs_B: (mitochondrial buffers equations)
-Eqs_CSQ: (calsequestrin equations)
-Eqs_NCX: (Na+/Ca2+ exchanger equations)
-Eqs_SOCE: (store operated Ca2+ entry equations)
-Eqs_NCE: (MITO Na+/Ca2+ exchanger equations)
-Eqs_MCU: (MITO Ca2+ uniporter equations)
-Eqs_MITO: (MITO equations (Eqs_B,Eqs_NCE and Eqs_MCU))
-Eqs_Diff: (Diffusion equations)
#-----------------------------------------------------------------------------
 Program 4. SCMI.py
#-----------------------------------------------------------------------------
-Description:generate the figures and data obtained with the single 
compartment model in fiber type I.
  
-Outputs:-outputs\caatps1.pdf (figure of [CaATP] during single AP in fiber
          type I with the SCM)
         -outputs\catns1.pdf (figure of [CaTn] during single AP in fiber
          type I with the SCM)
	 -outputs\cadyes1.pdf (figure of [CaDye] during single AP in fiber
          type I with the SCM)
	 -outputs\capvs1.pdf (figure of [CaPv] during single AP in fiber
          type I with the SCM)
	 -outputs\casercas1.pdf (figure of [CaSERCA] during single AP in 
	  fiber type I with the SCM)
	 -outputs\cancxs1.pdf (figure of [CaNCX] during single AP in fiber
          type I with the SCM)
	 -outputs\camitos1.pdf (figure of [Ca2+]mito during single AP in 
          fiber type I with the SCM)
	 -outputs\cats1.pdf (figure of [Ca2+]total during single AP in 
          fiber type I with the SCM)
	 -outputs\jrel1.pdf (figure of release rate of Ca2+ during single
          AP in fiber type I with the SCM)
	 -outputs\cs1.xls (table with the data that generate previus 
          figures)  
#-----------------------------------------------------------------------------
 Program 5. SCMIIA.py
#-----------------------------------------------------------------------------
-Description:generate the figures and data obtained with the single 
compartment model in fiber type IIA.
  
-Outputs:-outputs\caatps2a.pdf (figure of [CaATP] during single AP in fiber
          type IIA with the SCM)
         -outputs\catns2a.pdf (figure of [CaTn] during single AP in fiber
          type IIA with the SCM)
	 -outputs\cadyes2a.pdf (figure of [CaDye] during single AP in fiber
          type IIA with the SCM)
	 -outputs\capvs2a.pdf (figure of [CaPv] during single AP in fiber
          type IIA with the SCM)
	 -outputs\casercas2a.pdf (figure of [CaSERCA] during single AP in 
	  fiber type IIA with the SCM)
	 -outputs\cancxs2a.pdf (figure of [CaNCX] during single AP in fiber
          type IIA with the SCM)
	 -outputs\camitos2a.pdf (figure of [Ca2+]mito during single AP in 
          fiber type IIA with the SCM)
	 -outputs\cats2a.pdf (figure of [Ca2+]total during single AP in 
          fiber type IIA with the SCM)
	 -outputs\jrel2a.pdf (figure of release rate of Ca2+ during single
          AP in fiber type IIA with the SCM) 
	 -outputs\cs2a.xls (table with the data that generate previus 
          figures) 
#-----------------------------------------------------------------------------
 Program 6. SCMIIX.py
#-----------------------------------------------------------------------------
-Description:generate the figures and data obtained with the single 
compartment model in fiber type IIX.
  
-Outputs:-outputs\caatps2x.pdf (figure of [CaATP] during single AP in fiber
          type IIX with the SCM)
         -outputs\catns2x.pdf (figure of [CaTn] during single AP in fiber
          type IIX with the SCM)
	 -outputs\cadyes2x.pdf (figure of [CaDye] during single AP in fiber
          type IIX with the SCM)
	 -outputs\capvs2x.pdf (figure of [CaPv] during single AP in fiber
          type IIX with the SCM)
	 -outputs\casercas2x.pdf (figure of [CaSERCA] during single AP in 
	  fiber type IIX with the SCM)
	 -outputs\cancxs2x.pdf (figure of [CaNCX] during single AP in fiber
          type IIX with the SCM)
	 -outputs\camitos2x.pdf (figure of [Ca2+]mito during single AP in 
          fiber type IIX with the SCM)
	 -outputs\cats2x.pdf (figure of [Ca2+]total during single AP in 
          fiber type IIX with the SCM)
	 -outputs\jrel2x.pdf (figure of release rate of Ca2+ during single
          AP in fiber type IIX with the SCM) 
	 -outputs\cs2x.xls (table with the data that generate previus 
          figures) 
#-----------------------------------------------------------------------------
 Program 7. SCMIIB.py
#-----------------------------------------------------------------------------
-Description:generate the figures and data obtained with the single 
compartment model in fiber type IIB.
  
-Outputs:-outputs\caatps2b.pdf (figure of [CaATP] during single AP in fiber
          type IIB with the SCM)
         -outputs\catns2b.pdf (figure of [CaTn] during single AP in fiber
          type IIB with the SCM)
	 -outputs\cadyes2b.pdf (figure of [CaDye] during single AP in fiber
          type IIB with the SCM)
	 -outputs\capvs2b.pdf (figure of [CaPv] during single AP in fiber
          type IIB with the SCM)
	 -outputs\casercas2b.pdf (figure of [CaSERCA] during single AP in 
	  fiber type IIB with the SCM)
	 -outputs\cancxs2b.pdf (figure of [CaNCX] during single AP in fiber
          type IIB with the SCM)
	 -outputs\camitos2b.pdf (figure of [Ca2+]mito during single AP in 
          fiber type IIB with the SCM)
	 -outputs\cats2b.pdf (figure of [Ca2+]total during single AP in 
          fiber type IIB with the SCM)
	 -outputs\jrel2b.pdf (figure of release rate of Ca2+ during single
          AP in fiber type IIB with the SCM) 
	 -outputs\cs2b.xls (table with the data that generate previus 
          figures) 
#-----------------------------------------------------------------------------
 Program 7. MCMIIB.py
#-----------------------------------------------------------------------------
-Description:generate the figures and data obtained with the multi
compartment model in fiber type IIB.
  
-Outputs:-outputs\caatpm2b.pdf (figure of [CaATP] during single AP in fiber
          type IIB with the MCM)
         -outputs\catnm2b.pdf (figure of [CaTn] during single AP in fiber
          type IIB with the MCM)
	 -outputs\cadyem2b.pdf (figure of [CaDye] during single AP in fiber
          type IIB with the MCM)
	 -outputs\capvm2b.pdf (figure of [CaPv] during single AP in fiber
          type IIB with the MCM)
	 -outputs\casercam2b.pdf (figure of [CaSERCA] during single AP in 
	  fiber type IIB with the MCM)
	 -outputs\cancxm2b.pdf (figure of [CaNCX] during single AP in fiber
          type IIB with the MCM)
	 -outputs\camitom2b.pdf (figure of [Ca2+]mito during single AP in 
          fiber type IIB with the MCM)
	 -outputs\catm2b.pdf (figure of [Ca2+]total during single AP in 
          fiber type IIB with the MCM)
	 -outputs\jrelm2b.pdf (figure of release rate of Ca2+ during single
          AP in fiber type IIB with the MCM) 
	 -outputs\cm2b.xls (table with the data (concentrations as a function
           of time) that generate previus figures)
#-----------------------------------------------------------------------------
 Program 7. MCMIIB.py
#-----------------------------------------------------------------------------
-Description:generate the figures and data obtained with the multi
compartment model in fiber type IIB.
  
-Outputs:-outputs\caatpm2b.pdf (figure of [CaATP] during single AP in fiber
          type IIB with the MCM)
         -outputs\catnm2b.pdf (figure of [CaTn] during single AP in fiber
          type IIB with the MCM)
	 -outputs\cadyem2b.pdf (figure of [CaDye] during single AP in fiber
          type IIB with the MCM)
	 -outputs\capvm2b.pdf (figure of [CaPv] during single AP in fiber
          type IIB with the MCM)
	 -outputs\casercam2b.pdf (figure of [CaSERCA] during single AP in 
	  fiber type IIB with the MCM)
	 -outputs\cancxm2b.pdf (figure of [CaNCX] during single AP in fiber
          type IIB with the MCM)
	 -outputs\camitom2b.pdf (figure of [Ca2+]mito during single AP in 
          fiber type IIB with the MCM)
	 -outputs\catm2b.pdf (figure of [Ca2+]total during single AP in 
          fiber type IIB with the MCM)
	 -outputs\jrelm2b.pdf (figure of release rate of Ca2+ during single
          AP in fiber type IIB with the MCM) 
	 -outputs\cm2b.xls (table with the data (concentrations as a function
           of time) that generate previus figures)
#-----------------------------------------------------------------------------
 Program 7. MCMIIB.py
#-----------------------------------------------------------------------------
-Description:generate the figures and data obtained with the multi
compartment model in fiber type IIB.
  
-Outputs:-outputs\caatpm2b.pdf (figure of [CaATP] during single AP in fiber
          type IIB with the MCM)
         -outputs\catnm2b.pdf (figure of [CaTn] during single AP in fiber
          type IIB with the MCM)
	 -outputs\cadyem2b.pdf (figure of [CaDye] during single AP in fiber
          type IIB with the MCM)
	 -outputs\capvm2b.pdf (figure of [CaPv] during single AP in fiber
          type IIB with the MCM)
	 -outputs\casercam2b.pdf (figure of [CaSERCA] during single AP in 
	  fiber type IIB with the MCM)
	 -outputs\cancxm2b.pdf (figure of [CaNCX] during single AP in fiber
          type IIB with the MCM)
	 -outputs\camitom2b.pdf (figure of [Ca2+]mito during single AP in 
          fiber type IIB with the MCM)
	 -outputs\catm2b.pdf (figure of [Ca2+]total during single AP in 
          fiber type IIB with the MCM)
	 -outputs\jrelm2b.pdf (figure of release rate of Ca2+ during single
          AP in fiber type IIB with the MCM) 
	 -outputs\cm2b.xls (table with the data (concentrations as a function
           of time) that generate previus figures)
#-----------------------------------------------------------------------------
 Program 7. MCMIIB.py
#-----------------------------------------------------------------------------
-Description:generate the figures and data obtained with the multi
compartment model in fiber type IIB.
  
-Outputs:-outputs\caatpm2b.pdf (figure of [CaATP] during single AP in fiber
          type IIB with the MCM)
         -outputs\catnm2b.pdf (figure of [CaTn] during single AP in fiber
          type IIB with the MCM)
	 -outputs\cadyem2b.pdf (figure of [CaDye] during single AP in fiber
          type IIB with the MCM)
	 -outputs\capvm2b.pdf (figure of [CaPv] during single AP in fiber
          type IIB with the MCM)
	 -outputs\casercam2b.pdf (figure of [CaSERCA] during single AP in 
	  fiber type IIB with the MCM)
	 -outputs\cancxm2b.pdf (figure of [CaNCX] during single AP in fiber
          type IIB with the MCM)
	 -outputs\camitom2b.pdf (figure of [Ca2+]mito during single AP in 
          fiber type IIB with the MCM)
	 -outputs\catm2b.pdf (figure of [Ca2+]total during single AP in 
          fiber type IIB with the MCM)
	 -outputs\jrelm2b.pdf (figure of release rate of Ca2+ during single
          AP in fiber type IIB with the MCM) 
	 -outputs\cm2b.xls (table with the data (concentrations as a function
           of time) that generate previus figures)


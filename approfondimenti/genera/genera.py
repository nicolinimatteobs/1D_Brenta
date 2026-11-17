# Genera i PDF di approfondimento del Resoconto_modello_1D.md in approfondimenti/.
# Uso, dalla cartella del progetto:  python approfondimenti/genera/genera.py [numeri, es. 1 4 17]
# Serve reportlab (python -m pip install reportlab) e matplotlib.
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stile import Documento
import pdf_01_07, pdf_08_14, pdf_15_21

TUTTI = pdf_01_07.PDF + pdf_08_14.PDF + pdf_15_21.PDF

if __name__ == "__main__":
    scelti = {int(a) for a in sys.argv[1:]}
    for numero, file, titolo, sotto, sez, funz in TUTTI:
        if scelti and numero not in scelti: continue
        t0 = time.time()
        D = Documento(numero, file, titolo, sotto, sez)
        funz(D)
        D.salva()
        print(f"{file:45s} {time.time()-t0:5.1f} s")

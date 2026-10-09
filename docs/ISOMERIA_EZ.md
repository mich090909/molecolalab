# Isomeria geometrica E/Z — implementazione didattica v0.5

## Algoritmo e limiti

Il modulo `src/stereo.js` implementa `analyzeEZ`, `summarizeEZ` e `invertEZ`, funzioni indipendenti dal motore Avogadro. L'assegnazione E/Z riguarda **solo C=C aciclici**, con ciascun carbonio collegato mediante due legami semplici a due sostituenti **espliciti**, inclusi gli idrogeni. Viene richiesta valenza locale 4 (un legame doppio e due singoli), assenza di isotopi/cariche formali e connettività senza anelli. Se due ligandi sono indistinguibili, E/Z non è definito; se il confronto CIP supera le regole implementate, si segnala *non risolto*.

1. Le priorità CIP delle due ramificazioni su ciascun carbonio sono confrontate per numero atomico, propagando i livelli successivi e i nodi duplicati semplici per i legami multipli nel ligando (come già avviene per R/S). La procedura **non è** un'implementazione completa delle CIP Sequence Rules IUPAC.
2. Il vettore che unisce i carboni costituisce l'asse del doppio legame. I vettori dai carboni verso i sostituenti di priorità maggiore vengono proiettati sul piano perpendicolare all'asse.
3. Il prodotto scalare normalizzato positivo fra le due proiezioni assegna **Z**, negativo assegna **E**, a condizione che il valore assoluto sia almeno **0,85** e che su ogni carbonio i due ligandi risultino sul lato opposto dell'asse. In caso di geometria quasi ortogonale, allineata o localmente incoerente non viene assegnato un descrittore.
4. L'operazione **Genera isomero opposto** taglia logicamente la connettività C=C, determina il frammento su un lato e lo ruota rigidamente di **180°** intorno all'asse; verifica quindi l'inversione E↔Z. Non cambia gli ordini dei legami né le distanze intraframmento. Questo è un metodo di *costruzione geometrica* dell'altro isomero, **non una simulazione cinetica di una torsione del doppio legame**.

Il comando per la molecola speculare usato per R/S **non** produce E↔Z. E e Z sono relazioni fra sostituenti su due carboni, preservate da una riflessione globale.

## Verifiche indipendenti

Il file `test/fixtures/ez_rdkit_cases.json` include **180 strutture 3D** e descrittori E/Z calcolati in precedenza con RDKit 2025.09.4, comprendendo alcheni semplici, alogenoalcheni, sostituenti ossigenati e nitrilici e confronti di catene. I test `test/stereo_ez.test.mjs` verificano, per ciascuna struttura, descrittore, inversione andata/ritorno, conservazione delle lunghezze dei legami e import/export JSON/CJSON. Si aggiungono casi limite: etene (non stereogenico), C=O, C≡C, assenza di H espliciti, isotopi, anelli e geometrie alterate.

Questi controlli dimostrano concordanza **nei casi del test**, non equivalenza al toolkit completo CIP di RDKit e non validano casi non inclusi (anelli, isotopi, descrittori condizionati, cumuleni, atomi chiralmente distinti in parità di numero atomico).

Le coordinate dei modelli di esempio sono state calcolate **offline** con RDKit (ETKDG/MMFF94). L'app distribuita non include RDKit, Avogadro, Qt o altri binari di calcolo chimico.

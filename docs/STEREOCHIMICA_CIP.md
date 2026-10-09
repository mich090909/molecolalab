# Guida a R/S — MolecolaLab Web 0.4

## Come eseguire una prima prova

Apri `standalone.html`, seleziona **(R)-butan-2-olo**, poi clicca sul centro C2. Il pannello *Configurazione assoluta · R/S* mostra `C2: R`, la sequenza di priorità `1: O`, `2: C (etile)`, `3: C (metile)`, `4: H`. I colori e le dimensioni delle sfere sono **rappresentazioni grafiche** e non modificano le priorità CIP.

Clicca **Crea la molecola speculare**: la configurazione diventa **S**. Se selezioni **(2R,3R)-acido tartarico**, la riflessione trasforma entrambi i centri in **S**. Il pulsante non inverte un solo centro: rispecchia l'intera molecola. Per ripristinare l'originale premi **Annulla**.

## Cosa viene effettivamente calcolato

1. Si riconosce un atomo **C o Si neutro, oppure N⁺**, con **quattro legami singoli** verso quattro atomi esplicitamente rappresentati.
2. Le priorità 1–4 sono calcolate dalla **numerazione atomica**; quando i primi atomi sono uguali, vengono confrontati i livelli successivi delle catene. I doppi/tripli legami nelle ramificazioni sono rappresentati mediante nodi terminali duplicati, su **grafi aciclici**.
3. Se i quattro sostituenti sono distinguibili e le coordinate sono genuinamente 3D, il segno dell'orientazione spaziale rispetto alle quattro priorità determina **R o S**.
4. La semplice rotazione della **camera** non modifica il descrittore. Una riflessione delle coordinate nello spazio, invece, inverte la configurazione. Le misure geometriche (distanze e angoli) non cambiano nella riflessione.

**Attenzione:** questo è un *sottoinsieme didattico* delle regole di Cahn–Ingold–Prelog, **non l'intero insieme di regole IUPAC**. Non usa un motore di chimica computazionale o una libreria di percezione di aromaticità.

## Situazioni deliberatamente non assegnate

- Anelli e percorsi CIP ciclici: la gestione corretta delle duplicazioni anulari richiede un digrafo CIP più completo.
- Isotopi (D, T, ¹³C ecc.) e masse atomiche isotopiche.
- H **impliciti**: completa o disegna esplicitamente gli H prima di provare un'attribuzione.
- Assi/elicità chirali, configurazioni E/Z, centri pseudoasimmetrici, regole CIP che dipendono da altre unità stereogeniche, aromaticità/risonanza non rappresentata correttamente.
- Geometrie coplanari o quasi degeneri; coordinate sovrapposte; valenze o formule chimiche non coerenti.
- Sistemi in cui le priorità rimangono indistinguibili con il confronto implementato.

Se un caso non è supportato, **nessun R/S viene inventato**. La presenza del badge R/S non certifica da sola l'accuratezza chimica dell'intero modello importato o costruito.

## Esempi verificati con RDKit

| Molecola | Centro/i e risultato atteso |
| --- | --- |
| (R)-bromoclorofluorometano | C2: **R** (Br > Cl > F > H) |
| (S)-bromoclorofluorometano | C2: **S** |
| (R)-butan-2-olo | C2: **R** |
| (S)-butan-2-olo | C2: **S** |
| (S)-acido lattico | C2: **S** (OH > COOH > CH₃ > H) |
| (S)-alanina | C2: **S** (NH₂ > COOH > CH₃ > H) |
| Acido tartarico (2R,3R) | C4: **R**, C6: **R** (ID dell'editor; posizioni 2 e 3 IUPAC) |

Coordinate generate *offline* con **RDKit 2025.09.4**, embedding ETKDG (seed 55), geometrie MMFF94 e stereodescrittori verificati con `Chem.AssignStereochemistryFrom3D`. L'applicazione **non distribuisce né esegue RDKit**.

Un ulteriore set di **300 molecole organiche acicliche** generato con la stessa metodologia è incluso in `test/fixtures/cip_rdkit_300.json`: confronti **300/300** sulle configurazioni di riferimento. Ciò **non equivale** a validazione universale del sottosistema CIP.

## Unità di misura e salvataggi

Le lunghezze sono calcolate in Ångström e gli angoli in gradi dalle coordinate correnti. Per archiviare una molecola con la propria geometria usa il progetto **JSON nativo** oppure `.cjson` (controllando i campi mancanti). L'assegnazione R/S viene **ricalcolata**, non letta da etichette stereo presenti in file MOL stereochimici: controlla il risultato dopo ogni importazione.


## Isomeria geometrica E/Z (novità 0.5)

Nel pannello E/Z le priorità CIP sono confrontate separatamente su ciascuno dei due carboni C=C con due ligandi espliciti. Z indica priorità maggiori dallo stesso lato; E indica priorità maggiori da lati opposti. La rotazione rigidamente costruita del frammento genera l'altro isomero, non simula una rotazione fisica intorno al doppio legame. Un'operazione speculare globale lascia E/Z invariato. Dettagli: [Isomeria E/Z](ISOMERIA_EZ.md) e [Guida completa](GUIDA_COMPLETA.md).
